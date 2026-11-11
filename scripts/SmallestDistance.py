#!/usr/bin/env python3
"""
Plots cumulative distribution functions of the smallest distances computed by
`SmallestDistance.java`:
- INS to a call of any type;
- INS to INS;
- DEL to a call of any type;
- DEL to DEL.

A single plot shows the curves for a TR BED and a non-TR BED of >=20bp calls,
and for a TR BED and a non-TR BED of >=50bp calls (INS in shades of red, DEL in
shades of blue; TR curves use a lighter shade than non-TR curves; >=20bp curves
are solid, >=50bp curves are dashed). Each curve gives, for every X bin, the
fraction of all calls of that type whose distance is <= X. Calls with no call of
the given type within the max distance (value -1) count in the denominator but
never in the numerator, so a curve may not reach 1; their number is reported in
the legend.

The X axis is limited to 0..MAX_X and bins are BIN_WIDTH bp long.

Usage: SmallestDistance.py output.png title tr_20.bed non_tr_20.bed tr_50.bed non_tr_50.bed

`title` is used verbatim as the title of the plot.

All inputs are BEDs produced by `SmallestDistance.java`, with format:

CHROM, first, last+1, ID, type, length,
to_any, to_any_overlaps, to_del, to_del_overlaps, to_ins, to_ins_overlaps,
to_dup, to_dup_overlaps, to_inv, to_inv_overlaps,
id_any, id_del, id_ins, id_dup, id_inv
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Same constants as in `SmallestDistance.java`
TYPE_DEL = 0
TYPE_INS = 1

# Histogram geometry
MAX_X = 1000
BIN_WIDTH = 10

# Columns of the BED
COL_TYPE = 4
COL_TO_ANY = 6
COL_TO_DEL = 8
COL_TO_INS = 10


# Colors: (INS, DEL)
COLORS_NON_TR = ("#b01010", "#1050b0")
COLORS_TR = ("#f4a0a0", "#9cc4ee")


def load(path):
    rows = []
    with open(path) as f:
        for line in f:
            tokens = line.rstrip("\n").split("\t")
            rows.append((int(tokens[COL_TYPE]), int(tokens[COL_TO_ANY]), int(tokens[COL_TO_DEL]), int(tokens[COL_TO_INS])))
    return np.array(rows, dtype=np.int64).reshape(-1, 4)


def plot_curve(ax, distances, label, color, linestyle):
    found = distances[distances >= 0]
    n_missing = len(distances) - len(found)
    edges = np.arange(0, MAX_X + BIN_WIDTH, BIN_WIDTH)
    counts, _ = np.histogram(found, bins=edges)
    # Fraction of ALL calls (including the ones with distance -1, which are never
    # <= X) whose distance falls in a bin ending at or before X.
    cdf = np.cumsum(counts) / max(len(distances), 1)
    ax.step(edges[1:], cdf, where="pre", color=color, linestyle=linestyle, linewidth=2, label=f"{label}")   #label=f"{label} (with neighbor={len(found)}, without neighbor={n_missing})")


def plot_file(ax, path, name, colors, linestyle):
    data = load(path)
    svtype, to_any, to_del, to_ins = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    is_ins = svtype == TYPE_INS
    is_del = svtype == TYPE_DEL
    #plot_curve(ax, to_any[is_ins], name + r", INS $\rightarrow$ any", colors[0], linestyle)
    plot_curve(ax, to_ins[is_ins], name + r", INS $\rightarrow$ INS", colors[0], linestyle)
    #plot_curve(ax, to_any[is_del], name + r", DEL $\rightarrow$ any", colors[1], linestyle)
    plot_curve(ax, to_del[is_del], name + r", DEL $\rightarrow$ DEL", colors[1], linestyle)


def main():
    output_png = sys.argv[1]
    title = sys.argv[2]
    tr_bed_20, non_tr_bed_20 = sys.argv[3], sys.argv[4]
    tr_bed_50, non_tr_bed_50 = sys.argv[5], sys.argv[6]

    fig, ax = plt.subplots(figsize=(5, 5))
    plot_file(ax, non_tr_bed_20, r"$\geq$20bp, non-TR", COLORS_NON_TR, "--")
    plot_file(ax, tr_bed_20, r"$\geq$20bp, TR", COLORS_TR, "--")
    plot_file(ax, non_tr_bed_50, r"$\geq$50bp, non-TR", COLORS_NON_TR, "-")
    plot_file(ax, tr_bed_50, r"$\geq$50bp, TR", COLORS_TR, "-")
    ax.set_title(title)
    ax.set_xlim(0, MAX_X)
    ax.set_ylim(0, 1)
    ax.set_box_aspect(1)  # Square axes
    ax.set_axisbelow(True)
    ax.grid(True, color="lightgray", linewidth=0.8)
    ax.set_xlabel("Smallest distance (bp)")
    ax.set_ylabel(r"Fraction of calls with distance $\leq$ d")
    ax.legend(loc="upper right", fontsize="small")
    fig.tight_layout()
    fig.savefig(output_png, dpi=150, bbox_inches="tight", pad_inches=0.05)


if __name__ == "__main__":
    main()
