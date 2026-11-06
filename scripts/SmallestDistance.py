#!/usr/bin/env python3
"""
Plots cumulative distribution functions of the smallest distances computed by
`SmallestDistance.java`:
- INS to a call of any type;
- INS to INS;
- DEL to a call of any type;
- DEL to DEL.

A single plot shows the four curves (INS in shades of red, DEL in shades of
blue; the lighter shade is "to any"). Each curve gives, for every X bin, the
fraction of all calls of that type whose distance is <= X. Calls with no call of
the given type within the max distance (value -1) count in the denominator but
never in the numerator, so a curve may not reach 1; their number is reported in
the legend.

The X axis is limited to 0..MAX_X and bins are BIN_WIDTH bp long.

Usage: plot_smallest_distance.py output.png PREFIX P Q R input1.bed [input2.bed ...]

PREFIX, P, Q, R are only used in the title of the plot: "PREFIX: Length>=P bp,
AF>=Q, POS distance <=R Mbp".

Every input is a BED produced by `SmallestDistance.java` (e.g. one per
chromosome), with format:

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


def load(paths):
    rows = []
    for path in paths:
        with open(path) as f:
            for line in f:
                tokens = line.rstrip("\n").split("\t")
                rows.append((int(tokens[COL_TYPE]), int(tokens[COL_TO_ANY]), int(tokens[COL_TO_DEL]), int(tokens[COL_TO_INS])))
    return np.array(rows, dtype=np.int64).reshape(-1, 4)


def plot_curve(ax, distances, label, color):
    found = distances[distances >= 0]
    n_missing = len(distances) - len(found)
    edges = np.arange(0, MAX_X + BIN_WIDTH, BIN_WIDTH)
    counts, _ = np.histogram(found, bins=edges)
    # Fraction of ALL calls (including the ones with distance -1, which are never
    # <= X) whose distance falls in a bin ending at or before X.
    cdf = np.cumsum(counts) / max(len(distances), 1)
    ax.step(edges[1:], cdf, where="pre", color=color, linewidth=2,
            label=f"{label} (with neighbor={len(found)}, without neighbor={n_missing})")


def main():
    output_png = sys.argv[1]
    title_prefix = sys.argv[2]
    min_length, min_af, max_pos_distance_mbp = sys.argv[3], sys.argv[4], sys.argv[5]
    data = load(sys.argv[6:])
    svtype, to_any, to_del, to_ins = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    is_ins = svtype == TYPE_INS
    is_del = svtype == TYPE_DEL

    fig, ax = plt.subplots(figsize=(8, 8))
    plot_curve(ax, to_any[is_ins], r"INS $\rightarrow$ any", "#f4a0a0")
    plot_curve(ax, to_ins[is_ins], r"INS $\rightarrow$ INS", "#b01010")
    plot_curve(ax, to_any[is_del], r"DEL $\rightarrow$ any", "#9cc4ee")
    plot_curve(ax, to_del[is_del], r"DEL $\rightarrow$ DEL", "#1050b0")
    ax.set_title(r"$\geq$" + min_length + r"bp, AF$\geq$" + min_af + r", POS distance $\leq$" + max_pos_distance_mbp + "Mbp, " + title_prefix)
    ax.set_xlim(0, MAX_X)
    ax.set_ylim(0, 1)
    ax.set_axisbelow(True)
    ax.grid(True, color="lightgray", linewidth=0.8)
    ax.set_xlabel("Smallest distance (bp)")
    ax.set_ylabel(r"Fraction of calls with distance $\leq$ d")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(output_png, dpi=150)


if __name__ == "__main__":
    main()
