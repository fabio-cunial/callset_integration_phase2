#!/usr/bin/env python3
"""Build a (SVLEN bin) x (sample) matrix of non-ref call counts from a
multi-sample VCF, and save it as a heatmap image.

Cell (i,j) = number of records whose genotype in sample j contains at least one
ALT allele (het or hom-alt) and whose |SVLEN| falls in the i-th bin. Records
shorter than 10kb are ignored, since the smallest bin starts there.

The VCF is scanned once, top to bottom. The file is read as bytes and the
per-record sample block is never split into Python strings: length-filtering
happens first, and the genotypes of surviving records are scanned with numpy at
C speed. That keeps cohort-scale VCFs (>10k samples) practical.
"""

import sys
import gzip
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")          # write-to-file only; no display needed
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

TAB = 9
COLON = 58
DIGIT_1 = 49
DIGIT_9 = 57


def build_bins():
    """Return (edges, labels) for the fixed bin set.

    Bins are [10kb,20kb), ..., [90kb,100kb), [100kb,200kb), ..., [400kb,500kb),
    [500kb,1mbp), >=1mbp. `edges` has one more element than `labels`; the last
    edge is infinite.
    """
    edges = [10000 * i for i in range(1, 11)]          # 10kb ... 100kb
    edges += [100000 * i for i in range(2, 6)]         # 200kb ... 500kb
    edges += [1000000, float("inf")]                   # 1mbp, open-ended

    def fmt(x):
        if x >= 1000000:
            return "%gmbp" % (x / 1000000)
        return "%gkb" % (x / 1000)

    labels = ["[%s,%s)" % (fmt(edges[i]), fmt(edges[i + 1]))
              for i in range(len(edges) - 2)]
    labels.append(">=%s" % fmt(edges[-2]))
    return edges, labels


def open_vcf(path):
    """Open a plain or bgzipped/gzipped VCF in binary mode; '-' means stdin."""
    if path == "-":
        return sys.stdin.buffer
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def get_svlen(info, ref, alt):
    """Return |SVLEN| for a record (all args are bytes), or None if unknown.

    Uses INFO/SVLEN when present; otherwise falls back to the REF/ALT length
    difference, which is only meaningful for sequence-resolved calls (symbolic
    ALTs and breakends without SVLEN are skipped).
    """
    start = info.find(b"SVLEN=")
    while start >= 0:
        if start == 0 or info[start - 1] == ord(";"):   # a real INFO key
            end = info.find(b";", start)
            value = info[start + 6:end if end >= 0 else len(info)]
            value = value.split(b",")[0]
            try:
                return abs(int(value))
            except ValueError:
                try:
                    return abs(int(float(value)))
                except ValueError:
                    return None
        start = info.find(b"SVLEN=", start + 1)
    if alt.startswith(b"<") or b"[" in alt or b"]" in alt:
        return None
    return abs(len(alt) - len(ref))


def alt_mask(sample_block, n_samples):
    """Return a boolean array: True where that sample's GT has a non-ref allele.

    `sample_block` is the raw bytes of the record from the first sample column
    to the end of the line (trailing newline included, harmlessly). A sample is
    ALT iff its GT sub-field -- everything before the first ':' of the column,
    per the VCF spec, which requires GT to come first when present -- contains a
    digit 1-9. That covers '0/1', '1|1', '1/2', multi-digit allele indices and
    haploid calls, while '0/0', './.' and '.' stay false.

    The whole block is scanned with vector ops instead of one Python call per
    sample, which is what makes the 12k-sample case tractable.
    """
    array = np.frombuffer(sample_block, dtype=np.uint8)
    is_tab = array == TAB
    starts = np.empty(int(is_tab.sum()) + 1, dtype=np.int64)
    starts[0] = 0
    np.add(np.flatnonzero(is_tab), 1, out=starts[1:])
    if len(starts) != n_samples:
        raise RuntimeError("Record has %d sample columns, header declares %d."
                           % (len(starts), n_samples))

    # A character belongs to its column's GT sub-field iff no ':' has appeared
    # since the column started, i.e. the running colon count still equals the
    # count at the column's first character.
    colons = np.cumsum(array == COLON, dtype=np.int32)
    at_start = np.zeros(len(starts), dtype=np.int32)
    at_start[1:] = colons[starts[1:] - 1]
    column_of = np.cumsum(is_tab, dtype=np.int32)
    in_gt = colons == at_start[column_of]

    non_ref = (array >= DIGIT_1) & (array <= DIGIT_9) & in_gt
    return np.bitwise_or.reduceat(non_ref, starts)


def build_matrix(vcf_path, progress_every=0):
    """Stream the VCF once and return (matrix, sample_names, bin_labels)."""
    edges, labels = build_bins()
    n_bins = len(labels)
    min_length = edges[0]
    samples = None
    matrix = None
    n_records = 0

    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith(b"#"):
                if line.startswith(b"##"):
                    continue
                samples = line.rstrip().decode().split("\t")[9:]
                matrix = np.zeros((n_bins, len(samples)), dtype=np.int64)
                continue
            if matrix is None:
                raise RuntimeError("Malformed VCF: no #CHROM header line.")

            # Split off only the 9 leading columns; the sample block stays one
            # bytes object rather than becoming 12k Python strings.
            columns = line.split(b"\t", 9)
            if len(columns) < 10:
                continue
            svlen = get_svlen(columns[7], columns[3], columns[4])
            if svlen is None or svlen < min_length:
                continue
            row = np.searchsorted(edges, svlen, side="right") - 1
            if row >= n_bins:          # can only happen for the open-ended bin
                row = n_bins - 1
            matrix[row] += alt_mask(columns[9], len(samples))

            n_records += 1
            if progress_every and n_records % progress_every == 0:
                sys.stderr.write("%d records binned\n" % n_records)
                sys.stderr.flush()

    if matrix is None:
        raise RuntimeError("Malformed VCF: no #CHROM header line.")
    return matrix, samples, labels


def plot_matrix(matrix, samples, labels, log_scale, image_path, dpi):
    """Draw the matrix as a heatmap and write it to `image_path`.

    The format is taken from the file extension (.png, .pdf, .svg, ...).
    """
    n_samples = len(samples)
    # Grow with the sample count, but stay inside a window that fits a screen:
    # at cohort scale the columns become a dense band, which is the point.
    figure_width = min(20.0, max(6.0, 0.35 * n_samples + 3.0))
    figure, axes = plt.subplots(figsize=(figure_width, 6.0))

    norm = None
    if log_scale and matrix.max() > 0:
        # Counts decay by orders of magnitude across the bins; a log ramp keeps
        # the long bins from washing out. Zeros are masked to the "bad" color.
        norm = LogNorm(vmin=1, vmax=max(matrix.max(), 2))
    # A flame ramp: dark red -> red -> orange -> yellow -> white as the count
    # grows. "hot" starts just above black, so zeros are masked out and painted
    # pure black instead, keeping "no calls" unmistakably distinct from "few".
    flame = plt.get_cmap("hot").copy()
    flame.set_bad("black")
    axes.set_facecolor("black")
    image = axes.imshow(np.ma.masked_equal(matrix, 0), aspect="auto",
                        cmap=flame, norm=norm)

    if n_samples <= 60:
        axes.set_xticks(range(n_samples))
        axes.set_xticklabels(samples, rotation=90, fontsize=7)
    axes.set_yticks(range(len(labels)))
    axes.set_yticklabels(labels, fontsize=8)
    axes.set_xlabel("sample (n=%d)" % n_samples)
    axes.set_ylabel("SVLEN bin")
    axes.set_title("Non-ref calls per SVLEN bin per sample")

    axes.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    if n_samples <= 60:            # cell borders only while cells are visible
        axes.set_xticks(np.arange(-0.5, n_samples, 1), minor=True)
    axes.grid(which="minor", color="black", linewidth=1.0)
    axes.tick_params(which="minor", length=0)

    colorbar = figure.colorbar(image, ax=axes)
    colorbar.set_label("number of records" + (" (log scale)" if norm else ""))

    # With few samples the exact counts fit inside the cells.
    if n_samples <= 20:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix[i, j]
                if value == 0:
                    continue
                # The flame gets brighter with the count, so the ink has to go
                # the other way: black on the yellow/white end, white on the
                # dark-red end. Measured on the same scale used for the fill.
                bright = np.asarray(image.norm(value)).item() > 0.55
                axes.text(j, i, str(value), ha="center", va="center",
                          fontsize=7, color="black" if bright else "white")

    figure.tight_layout()
    figure.savefig(image_path, dpi=dpi)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="multi-sample VCF (plain or gzipped; - for stdin)")
    parser.add_argument("image", help="output heatmap file; the extension picks "
                                      "the format (.png, .pdf, .svg, ...)")
    parser.add_argument("--tsv", help="also write the matrix to this TSV file")
    parser.add_argument("--linear", action="store_true",
                        help="use a linear color scale instead of logarithmic")
    parser.add_argument("--dpi", type=int, default=150,
                        help="resolution of the output image (default: 150)")
    parser.add_argument("--progress", type=int, default=0, metavar="N",
                        help="report to stderr every N binned records")
    args = parser.parse_args()

    matrix, samples, labels = build_matrix(args.vcf, args.progress)
    if args.tsv:
        with open(args.tsv, "wt") as out:
            out.write("bin\t" + "\t".join(samples) + "\n")
            for label, row in zip(labels, matrix):
                out.write(label + "\t" + "\t".join(str(x) for x in row) + "\n")
    plot_matrix(matrix, samples, labels, not args.linear, args.image, args.dpi)


if __name__ == "__main__":
    main()
