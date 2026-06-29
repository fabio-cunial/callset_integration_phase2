#!/usr/bin/env python3
"""
WARNING: this script is AI-generated.

Plot an assembly-to-reference alignment TSV as a "matrix of arrows".

Input TSV (tab-separated, no header), one row per alignment, as produced by
`StudyAssemblySam.java`:

    contigID  contigLength  alignmentStartInContig  alignmentEndInContig
    alignedChromosomeId  alignmentStartInChromosome  alignmentEndInChromosome
    forwardOrReverse(FW|RC)  MAPQ

Each contig gets its own horizontal row. Every alignment of that contig is
drawn as an arrow whose horizontal extent is [alignmentStartInContig,
alignmentEndInContig], whose direction encodes strand (FW -> right, RC ->
left), and whose color encodes the aligned reference chromosome. A legend maps
colors to chromosomes.

Usage:
    python3 plot_assembly_arrows.py input.tsv output.png
    python3 plot_assembly_arrows.py input.tsv output.png --min-mapq 20
"""

import argparse
import csv
import math
import os
import sys
from collections import OrderedDict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
from matplotlib.lines import Line2D


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_tsv", help="Input TSV file (StudyAssemblySam output).")
    ap.add_argument("output_image", help="Output image file (.png/.pdf/.svg).")
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="Skip alignments with MAPQ below this value (default: 0).")
    ap.add_argument("--row-height", type=float, default=0.6,
                    help="Vertical space allotted to each contig row (default: 0.6).")
    ap.add_argument("--width", type=float, default=14.0,
                    help="Figure width in inches (default: 14).")
    ap.add_argument("--dpi", type=int, default=150,
                    help="Output resolution (default: 150).")
    ap.add_argument("--sort", choices=["input", "length", "name"], default="input",
                    help="Order of contig rows top-to-bottom (default: input order).")
    ap.add_argument("--hist-tsv", default=None,
                    help="Optional TSV (tab-separated, no header) with columns "
                         "'length<TAB>type<TAB>mapq<TAB>...', where type is an "
                         "integer label. If given, an extra histogram plot is written.")
    ap.add_argument("--hist-output", default=None,
                    help="Output image for the per-type length histogram "
                         "(required if --hist-tsv is given).")
    ap.add_argument("--hist-mapq-output", default=None,
                    help="Output image for the per-type MAPQ histogram. If "
                         "omitted, it is derived from --hist-output by inserting "
                         "'_mapq' before the extension.")
    ap.add_argument("--hist-bins", type=int, default=50,
                    help="Number of histogram bins (default: 50).")
    return ap.parse_args()


# Label used for any chromosome that is not one of the canonical human ones.
OTHER_LABEL = "other"


# Human-readable names for the integer 'type' labels in the --hist-tsv input.
CATEGORY_NAMES = {
    0: "most frequent chr and orientation, concordant pos",
    1: "most frequent chr and orientation, discordant pos",
    2: "most frequent chr, opposite orientation",
    3: "different chr",
}


def canonical_rank(chrom):
    """Return the canonical rank of a chromosome, or None if non-canonical.

    Canonical = chr1..chr22, chrX, chrY, chrM (with or without a 'chr' prefix).
    """
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        n = int(name)
        return n if 1 <= n <= 22 else None
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return special.get(name.upper())


def color_key(chrom):
    """The chromosome's color bucket: itself if canonical, else OTHER_LABEL."""
    return chrom if canonical_rank(chrom) is not None else OTHER_LABEL


def chromosome_sort_key(label):
    """Order color buckets: chr1..chr22, chrX, chrY, chrM, then 'other' last."""
    if label == OTHER_LABEL:
        return (2, 0, "")
    rank = canonical_rank(label)
    if rank is not None:
        return (0, rank, "")
    return (1, 0, label)


def load_alignments(path, min_mapq):
    """Return (contigs, rows).

    contigs: OrderedDict contigID -> contigLength (first-seen order)
    rows:    list of dicts with alignment fields.
    """
    contigs = OrderedDict()
    rows = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for lineno, fields in enumerate(reader, start=1):
            if not fields or (len(fields) == 1 and not fields[0].strip()):
                continue
            if len(fields) < 9:
                sys.stderr.write(
                    f"Warning: line {lineno} has {len(fields)} columns (<9), skipping.\n")
                continue
            try:
                contig_id = fields[0]
                contig_len = int(fields[1])
                start_in_contig = int(fields[2])
                end_in_contig = int(fields[3])
                chrom = fields[4]
                strand = fields[7].strip().upper()
                mapq = int(fields[8])
            except ValueError as e:
                sys.stderr.write(f"Warning: line {lineno} failed to parse ({e}), skipping.\n")
                continue
            if mapq < min_mapq:
                continue
            contigs.setdefault(contig_id, contig_len)
            # contigLength should be consistent across a contig's alignments;
            # keep the max just in case of clipping discrepancies.
            if contig_len > contigs[contig_id]:
                contigs[contig_id] = contig_len
            rows.append({
                "contig": contig_id,
                "start": min(start_in_contig, end_in_contig),
                "end": max(start_in_contig, end_in_contig),
                "chrom": chrom,
                "reverse": strand == "RC",
                "mapq": mapq,
            })
    return contigs, rows


def build_color_map(rows):
    """Assign a color to each chromosome bucket, in natural order.

    Every canonical chromosome (chr1..chr22, chrX, chrY, chrM) gets its own
    color; all non-canonical chromosomes share a single 'other' color/legend
    entry.
    """
    buckets = sorted({color_key(r["chrom"]) for r in rows}, key=chromosome_sort_key)
    canonical = [b for b in buckets if b != OTHER_LABEL]
    # tab20 gives 20 distinguishable colors; fall back to hsv for more.
    if len(canonical) <= 20:
        cmap = plt.get_cmap("tab20")
        colors = [cmap(i % 20) for i in range(len(canonical))]
    else:
        cmap = plt.get_cmap("hsv")
        n = max(len(canonical), 1)
        colors = [cmap(i / n) for i in range(len(canonical))]
    color_map = OrderedDict(zip(canonical, colors))
    if OTHER_LABEL in buckets:
        color_map[OTHER_LABEL] = (0.5, 0.5, 0.5, 1.0)  # gray for everything else
    return color_map


def order_contigs(contigs, how):
    items = list(contigs.items())
    if how == "length":
        items.sort(key=lambda kv: kv[1], reverse=True)
    elif how == "name":
        items.sort(key=lambda kv: kv[0])
    # "input" keeps first-seen order
    return [c for c, _ in items]


def load_metrics_by_type(path):
    """Return an OrderedDict type(int) -> {"length": [float], "mapq": [float]}.

    Reads a tab-separated file with columns 'length<TAB>type<TAB>mapq<TAB>...';
    any further columns are ignored. Types are kept in first-seen order.
    """
    by_type = OrderedDict()
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for lineno, fields in enumerate(reader, start=1):
            if not fields or (len(fields) == 1 and not fields[0].strip()):
                continue
            if len(fields) < 3:
                sys.stderr.write(
                    f"Warning: hist line {lineno} has {len(fields)} columns (<3), skipping.\n")
                continue
            try:
                length = float(fields[0])
                type_label = int(fields[1])
                mapq = float(fields[2])
            except ValueError as e:
                sys.stderr.write(
                    f"Warning: hist line {lineno} failed to parse ({e}), skipping.\n")
                continue
            metrics = by_type.setdefault(type_label, {"length": [], "mapq": []})
            metrics["length"].append(length)
            metrics["mapq"].append(mapq)
    return by_type


def plot_histograms_by_type(values_by_type, output_image, n_bins, width, dpi,
                            xlabel, log_x, log_y=False):
    """Plot the histogram of one metric for each type in its own panel.

    `values_by_type` maps type(int) -> list of values. All panels share the
    same X axis and bin edges so the distributions are directly comparable
    across categories. With `log_x=True` the axis is log-scaled and the bins
    are logarithmically spaced (non-positive values are dropped); otherwise the
    axis and bins are linear. With `log_y=True` the count (Y) axis is
    log-scaled.
    """
    all_vals = [v for vals in values_by_type.values() for v in vals]
    if log_x:
        # Positive values only, with logarithmically-spaced bin edges so bars
        # are uniform-width on the log axis.
        usable = [v for v in all_vals if v > 0]
        n_dropped = len(all_vals) - len(usable)
        if n_dropped:
            sys.stderr.write(
                f"Warning: dropped {n_dropped} non-positive {xlabel} value(s) "
                f"(incompatible with a log X axis).\n")
    else:
        usable = list(all_vals)
    if not usable:
        sys.stderr.write(f"No {xlabel} values to plot. Nothing written to {output_image}.\n")
        return
    # Shared bin edges so the per-panel histograms are directly comparable.
    lo, hi = min(usable), max(usable)
    if log_x:
        if lo == hi:
            hi = lo * 10.0  # avoid a degenerate single-edge range
        log_lo, log_hi = math.log10(lo), math.log10(hi)
        bins = [10.0 ** (log_lo + (log_hi - log_lo) * k / n_bins)
                for k in range(n_bins + 1)]
    else:
        if lo == hi:
            hi = lo + 1.0  # avoid a degenerate single-edge range
        bins = [lo + (hi - lo) * k / n_bins for k in range(n_bins + 1)]

    types = sorted(values_by_type)
    n_panels = len(types)
    fig, axes = plt.subplots(n_panels, 1, sharex=True,
                             figsize=(width, max(2.0, 2.2 * n_panels)),
                             squeeze=False)
    for ax, type_label in zip(axes[:, 0], types):
        vals = [v for v in values_by_type[type_label] if (v > 0 or not log_x)]
        ax.hist(vals, bins=bins, color="C0")
        if log_x:
            ax.set_xscale("log")
        if log_y:
            ax.set_yscale("log")
        ax.set_ylabel("count")
        title = CATEGORY_NAMES.get(type_label, f"type {type_label}")
        ax.set_title(f"{title}  (n={len(vals)})", fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[-1, 0].set_xlabel(xlabel)

    fig.tight_layout()
    fig.savefig(output_image, dpi=dpi, bbox_inches="tight")
    sys.stderr.write(
        f"Wrote {output_image}: {n_panels} panels, {len(all_vals)} values.\n")


def main():
    args = parse_args()

    if args.hist_tsv:
        if not args.hist_output:
            sys.stderr.write("Error: --hist-output is required when --hist-tsv is given.\n")
            sys.exit(1)
        by_type = load_metrics_by_type(args.hist_tsv)
        if not by_type:
            sys.stderr.write("No histogram values to plot. Nothing written.\n")
        else:
            length_by_type = OrderedDict(
                (t, m["length"]) for t, m in by_type.items())
            plot_histograms_by_type(length_by_type, args.hist_output,
                                    args.hist_bins, args.width, args.dpi,
                                    xlabel="length", log_x=True, log_y=True)
            mapq_output = args.hist_mapq_output
            if mapq_output is None:
                root, ext = os.path.splitext(args.hist_output)
                mapq_output = f"{root}_mapq{ext}"
            mapq_by_type = OrderedDict(
                (t, m["mapq"]) for t, m in by_type.items())
            plot_histograms_by_type(mapq_by_type, mapq_output,
                                    args.hist_bins, args.width, args.dpi,
                                    xlabel="MAPQ", log_x=False, log_y=True)
    contigs, rows = load_alignments(args.input_tsv, args.min_mapq)

    if not rows:
        sys.stderr.write("No alignments to plot (after filtering). Nothing written.\n")
        sys.exit(1)

    color_map = build_color_map(rows)
    contig_order = order_contigs(contigs, args.sort)
    contig_y = {c: i for i, c in enumerate(contig_order)}

    # Group alignments by contig for drawing.
    by_contig = {c: [] for c in contig_order}
    for r in rows:
        by_contig[r["contig"]].append(r)

    n_contigs = len(contig_order)

    fig_height = max(2.0, n_contigs * args.row_height + 1.5)
    fig, ax = plt.subplots(figsize=(args.width, fig_height))

    arrow_thickness = args.row_height * 0.45
    # Every contig occupies the same horizontal space [0, 1]; positions are
    # expressed as a fraction of each contig's own length. Arrowhead length is
    # therefore a constant fraction of that shared width.
    head_len = 0.006

    for contig in contig_order:
        y = contig_y[contig]
        length = contigs[contig] or 1  # guard against zero-length contigs
        # Faint baseline spanning the whole (normalized) contig.
        ax.plot([0, 1], [y, y], color="0.85", linewidth=1.0, zorder=1)
        for aln in by_contig[contig]:
            start = aln["start"] / length
            end = aln["end"] / length
            span = end - start
            if aln["reverse"]:
                x_tail, dx = end, -span
            else:
                x_tail, dx = start, span
            # Cap arrowhead length so very short alignments still render.
            this_head = min(head_len, abs(dx)) if dx != 0 else head_len
            color = color_map[color_key(aln["chrom"])]
            arrow = FancyArrow(
                x_tail, y, dx, 0,
                width=arrow_thickness * 0.5,
                head_width=arrow_thickness,
                head_length=this_head,
                length_includes_head=True,
                color=color,
                alpha=0.9,
                zorder=2,
            )
            ax.add_patch(arrow)

    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-1, n_contigs)
    ax.invert_yaxis()  # first contig at top
    ax.set_yticks(range(n_contigs))
    ax.set_yticklabels(contig_order, fontsize=8)
    ax.set_xlabel("Position in contig (fraction of contig length)")
    ax.set_ylabel("Contig")
    ax.set_title("Assembly-to-reference alignments")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend: color -> chromosome, plus a note on arrow direction.
    legend_handles = [
        Line2D([0], [0], color=color, lw=6, label=chrom)
        for chrom, color in color_map.items()
    ]
    ax.legend(
        handles=legend_handles,
        title="Aligned chromosome\n(arrow points 5'->3' on contig;\nleft-pointing = reverse strand)",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        fontsize=8,
        title_fontsize=8,
        ncol=1 if len(color_map) <= 26 else 2,
    )

    fig.tight_layout()
    fig.savefig(args.output_image, dpi=args.dpi, bbox_inches="tight")
    sys.stderr.write(
        f"Wrote {args.output_image}: {n_contigs} contigs, {len(rows)} alignments, "
        f"{len(color_map)} chromosomes.\n")


if __name__ == "__main__":
    main()
