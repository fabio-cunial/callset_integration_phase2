#!/usr/bin/env python3
"""Count BNDs per sample in a multi-sample VCF, split by the kind of contig pair
they join, and save the categories as a jittered scatterplot.

One point per sample per category: its y is the number of BND records that are
non-ref (het or hom-alt) in that sample, its x is the category position plus a
small random offset, so that samples with equal counts stay distinguishable.

A second panel, to the right of the first one and drawn the same way, breaks the
same records down by standard chromosome instead: one column per chr1...chr22,
chrX, chrY, chrM, and a record counts for each of the (one or two) standard
contigs it touches. The two panels have independent y axes.

A contig is *standard* when it is one of chr1...chr22, chrX, chrY, chrM (names
without the `chr` prefix and the `MT` alias are accepted too); everything else
-- `chr*_*_random`, `chrUn_*`, `chr*_alt`, `chrEBV`, decoys -- is non-standard.
The categories **overlap**: each record is counted in every one it matches, so
the columns do not add up to the total number of BNDs.

  intra-chromosomal             both contigs the same
  inter-chromosomal             the two contigs differ
  chrM + standard               one contig chrM, the other a different standard
  chrEBV + standard             one contig chrEBV, the other standard
  standard + non-standard       exactly one contig non-standard
  non-standard + non-standard   both contigs non-standard, same one included

A chrEBV-to-chr1 breakend, for instance, lands in inter-chromosomal, in
chrEBV + standard and in standard + non-standard at once.

The VCF is scanned once, top to bottom. It is read as bytes and the per-record
sample block is never split into Python strings: genotypes are scanned with
numpy at C speed, which keeps cohort-scale VCFs (>10k samples) practical.
"""

import sys
import gzip
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")          # write-to-file only; no display needed
import matplotlib.pyplot as plt

TAB = 9
COLON = 58
DIGIT_1 = 49
DIGIT_9 = 57

INTRA = 0
INTER = 1
M_STD = 2
EBV_STD = 3
STD_NONSTD = 4
NONSTD_NONSTD = 5
N_CATEGORIES = 6

# Single-line names, used for the TSV header.
CATEGORIES = ["intra-chromosomal", "inter-chromosomal", "chrM + standard",
              "chrEBV + standard", "standard + non-standard",
              "non-standard + non-standard"]
# Same order, wrapped so that six ticks fit under the axis.
CATEGORY_LABELS = ["intra-\nchromosomal", "inter-\nchromosomal",
                   "chrM +\nstandard", "chrEBV +\nstandard",
                   "standard +\nnon-standard", "non-standard +\nnon-standard"]
COLORS = ["#2f6db3", "#c2650f", "#3f8f5b", "#a04ba0", "#b3312f", "#7a6a55"]

MITO = b"chrM"
EBV = b"chrEBV"
# Column order of the per-chromosome panel; also the definition of "standard".
CHROMOSOMES = ([b"chr%d" % i for i in range(1, 23)]
               + [b"chrX", b"chrY", MITO])
CHROMOSOME_INDEX = {name: i for i, name in enumerate(CHROMOSOMES)}
N_CHROMOSOMES = len(CHROMOSOMES)
# The `chr` prefix is dropped: 25 ticks only fit if the labels are short.
CHROMOSOME_LABELS = [name.decode()[3:] for name in CHROMOSOMES]
# Alternated column by column, which is all that is needed to tell the narrow
# neighbouring clouds apart.
CHROMOSOME_COLORS = ["#2f6db3", "#79a8d8"]
STANDARD = frozenset(CHROMOSOMES)


def canonical(contig):
    """Return `contig` in `chr`-prefixed form, with `chrMT` folded onto chrM.

    Callsets built against the no-prefix flavour of GRCh38 name the contigs
    `1`, `X`, `MT`; normalising here lets one contig table serve both.
    """
    name = contig if contig.startswith(b"chr") else b"chr" + contig
    return MITO if name == b"chrMT" else name


def categories_of(chrom, mate):
    """Return every category index a BND joining `chrom` and `mate` falls into.

    The categories overlap, so this is a list, never a single index. Intra vs
    inter always contributes exactly one entry; how standard the two contigs
    are contributes zero, one or two more.
    """
    first = canonical(chrom)
    second = canonical(mate)
    first_standard = first in STANDARD
    second_standard = second in STANDARD

    found = [INTRA if first == second else INTER]
    if not first_standard and not second_standard:
        found.append(NONSTD_NONSTD)
    elif not first_standard or not second_standard:
        found.append(STD_NONSTD)
        if (first if not first_standard else second) == EBV:
            found.append(EBV_STD)
    elif (first == MITO) != (second == MITO):
        # Both standard and exactly one is chrM, so the mate is "another"
        # standard chromosome; a chrM-to-chrM breakend is intra-chromosomal.
        found.append(M_STD)
    return found


def chromosomes_of(chrom, mate):
    """Return the column indices of the standard contigs a BND touches.

    Zero, one or two entries: non-standard contigs have no column, and an
    intra-chromosomal breakend touches its chromosome once, not twice.
    """
    found = []
    for contig in {canonical(chrom), canonical(mate)}:
        column = CHROMOSOME_INDEX.get(contig)
        if column is not None:
            found.append(column)
    return found


def open_vcf(path):
    """Open a plain or bgzipped/gzipped VCF in binary mode; '-' means stdin."""
    if path == "-":
        return sys.stdin.buffer
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def mate_chrom(alt, info):
    """Return the contig of a BND's mate (bytes), or None if not determinable.

    The breakend ALT syntax puts the mate locus between a pair of brackets --
    't[chr2:1000[', ']chr2:1000]t' and the other two orientations -- so the
    contig is what precedes the last ':' inside the brackets. Symbolic <BND>
    records without brackets fall back to INFO/CHR2; single breakends ('.t',
    't.') have no mate and yield None.
    """
    first = alt.split(b",")[0]
    open_bracket = -1
    close_bracket = -1
    for bracket in (b"[", b"]"):
        position = first.find(bracket)
        if position >= 0 and (open_bracket < 0 or position < open_bracket):
            open_bracket = position
            close_bracket = first.find(bracket, position + 1)
    if open_bracket >= 0 and close_bracket > open_bracket:
        locus = first[open_bracket + 1:close_bracket]
        cut = locus.rfind(b":")
        if cut > 0:
            return locus[:cut]

    start = info.find(b"CHR2=")
    while start >= 0:
        if start == 0 or info[start - 1] == ord(";"):   # a real INFO key
            end = info.find(b";", start)
            return info[start + 5:end if end >= 0 else len(info)]
        start = info.find(b"CHR2=", start + 1)
    return None


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


def build_counts(vcf_path, progress_every=0):
    """Stream the VCF once, returning (counts, chrom_counts, samples, n_skipped).

    `counts` is an N_CATEGORIES x n_samples int64 array, one row per entry of
    CATEGORIES; a record contributes to every category it matches, so the rows
    overlap. `chrom_counts` is an N_CHROMOSOMES x n_samples int64 array in
    CHROMOSOMES order, holding the same records broken down by the standard
    contigs they touch. `n_skipped` counts records whose mate contig could not
    be determined (single breakends, or symbolic ALTs without INFO/CHR2); those
    contribute to neither array.
    """
    samples = None
    counts = None
    chrom_counts = None
    n_records = 0
    n_skipped = 0

    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith(b"#"):
                if line.startswith(b"##"):
                    continue
                samples = line.rstrip().decode().split("\t")[9:]
                counts = np.zeros((N_CATEGORIES, len(samples)), dtype=np.int64)
                chrom_counts = np.zeros((N_CHROMOSOMES, len(samples)),
                                        dtype=np.int64)
                continue
            if counts is None:
                raise RuntimeError("Malformed VCF: no #CHROM header line.")

            # Split off only the 9 leading columns; the sample block stays one
            # bytes object rather than becoming one Python string per sample.
            columns = line.split(b"\t", 9)
            if len(columns) < 10:
                continue
            mate = mate_chrom(columns[4], columns[7])
            if mate is None:
                n_skipped += 1
                continue
            # One genotype scan per record, added to each category it matches.
            mask = alt_mask(columns[9], len(samples))
            for category in categories_of(columns[0], mate):
                counts[category] += mask
            for column in chromosomes_of(columns[0], mate):
                chrom_counts[column] += mask

            n_records += 1
            if progress_every and n_records % progress_every == 0:
                sys.stderr.write("%d records counted\n" % n_records)
                sys.stderr.flush()

    if counts is None:
        raise RuntimeError("Malformed VCF: no #CHROM header line.")
    return counts, chrom_counts, samples, n_skipped


def draw_clouds(axes, matrix, colors, generator, size, alpha):
    """Scatter one jittered cloud per row of `matrix`, and return the medians.

    Column `i` of the panel holds row `i` of the matrix, drawn at x=i plus a
    random offset per sample; a short horizontal rule marks the median, the one
    summary worth reading off a cloud directly.
    """
    n_samples = matrix.shape[1]
    medians = []
    for column in range(matrix.shape[0]):
        values = matrix[column]
        jitter = generator.uniform(-0.18, 0.18, size=n_samples)
        axes.scatter(column + jitter, values, s=size, alpha=alpha,
                     color=colors[column % len(colors)], linewidths=0.0,
                     zorder=2)
        median = float(np.median(values)) if n_samples else 0.0
        axes.hlines(median, column - 0.3, column + 0.3, color="#333333",
                    linewidth=2.0, zorder=3)
        medians.append(median)
    return medians


def style_panel(axes, n_columns, log_scale):
    """Apply the shared look: y grid behind the points, no top/right spines."""
    axes.set_xlim(-0.6, n_columns - 0.4)
    if log_scale:
        axes.set_yscale("log")
    else:
        axes.set_ylim(bottom=0)
    axes.grid(axis="y", color="#dddddd", linewidth=0.8, zorder=0)
    axes.set_axisbelow(True)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)


def plot_counts(counts, chrom_counts, samples, image_path, dpi, seed,
                log_scale):
    """Draw the category panel and the per-chromosome one, write `image_path`.

    Both panels are jittered point clouds of the same quantity -- non-ref BND
    records per sample -- binned differently on x. Their y axes are independent:
    inter-chromosomal alone outweighs any single chromosome by enough that a
    shared axis would flatten the right panel into a band.
    The format is taken from the file extension (.png, .pdf, .svg, ...).
    """
    n_samples = len(samples)
    generator = np.random.default_rng(seed)
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(17.0, 6.0), gridspec_kw={"width_ratios": [3.0, 4.0]})

    # Denser clouds need smaller, more transparent points to stay readable.
    if n_samples <= 100:
        size, alpha = 42.0, 0.85
    elif n_samples <= 2000:
        size, alpha = 16.0, 0.5
    else:
        size, alpha = 6.0, 0.25

    medians = draw_clouds(left, counts, COLORS, generator, size, alpha)
    # With six clouds side by side there is no room for a label next to the
    # median line, so it goes under the category name instead.
    left.set_xticks(list(range(N_CATEGORIES)))
    left.set_xticklabels(["%s\nmedian %.0f" % (label, median)
                          for label, median in zip(CATEGORY_LABELS, medians)],
                         fontsize=9)
    left.set_xlabel("BND category (categories overlap; a record can be in several)")
    left.set_ylabel("number of non-ref BND records per sample")
    style_panel(left, N_CATEGORIES, log_scale)

    # 25 columns leave room for the chromosome name and nothing else, so the
    # medians here are only the drawn rules.
    draw_clouds(right, chrom_counts, CHROMOSOME_COLORS, generator, size, alpha)
    right.set_xticks(list(range(N_CHROMOSOMES)))
    right.set_xticklabels(CHROMOSOME_LABELS, fontsize=8)
    right.set_xlabel("chromosome (a record counts once for each standard "
                     "contig it touches)")
    right.set_ylabel("number of non-ref BND records per sample")
    style_panel(right, N_CHROMOSOMES, log_scale)

    figure.suptitle("BND calls per sample (n=%d)" % n_samples)
    figure.tight_layout()
    figure.savefig(image_path, dpi=dpi)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="multi-sample BND-only VCF (plain or gzipped; - for stdin)")
    parser.add_argument("image", help="output plot file; the extension picks "
                                      "the format (.png, .pdf, .svg, ...)")
    parser.add_argument("--tsv", help="also write the per-sample counts to this TSV file")
    parser.add_argument("--log", action="store_true",
                        help="use a logarithmic y axis")
    parser.add_argument("--dpi", type=int, default=150,
                        help="resolution of the output image (default: 150)")
    parser.add_argument("--seed", type=int, default=42,
                        help="seed of the x jitter, for reproducible plots (default: 42)")
    parser.add_argument("--progress", type=int, default=0, metavar="N",
                        help="report to stderr every N counted records")
    args = parser.parse_args()

    counts, chrom_counts, samples, n_skipped = build_counts(args.vcf,
                                                            args.progress)
    if n_skipped:
        sys.stderr.write("Warning: %d records skipped (no mate contig in ALT "
                         "or INFO/CHR2).\n" % n_skipped)
    if args.tsv:
        with open(args.tsv, "wt") as out:
            out.write("sample\t" + "\t".join(CATEGORIES) + "\t"
                      + "\t".join(name.decode() for name in CHROMOSOMES) + "\n")
            for j, sample in enumerate(samples):
                out.write(sample + "\t"
                          + "\t".join(str(counts[c, j])
                                      for c in range(N_CATEGORIES)) + "\t"
                          + "\t".join(str(chrom_counts[c, j])
                                      for c in range(N_CHROMOSOMES)) + "\n")
    plot_counts(counts, chrom_counts, samples, args.image, args.dpi, args.seed,
                args.log)


if __name__ == "__main__":
    main()
