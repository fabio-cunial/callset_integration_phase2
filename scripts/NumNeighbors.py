#!/usr/bin/env python3
"""
Plots, for every N in a set of bins, the fraction of calls that have at least N
neighbors. A single plot shows the curves for INS and DEL, TR and non-TR, >=20bp
and >=50bp calls (INS in shades of red, DEL in shades of blue; TR curves use a
lighter shade than non-TR curves; >=20bp curves are dashed, >=50bp curves are
solid).

The values of N are 1, 2, 4, 8, 16, 20, and then increments of 10 up to
MAX_N, which is also the upper limit of the X axis. The X axis is in log scale
(from 1 to MAX_N) and the Y axis is linear.

Usage: NumNeighbors.py output.png title [indir]

`title` is used verbatim as the title of the plot. `indir` (default: the current
directory) must contain the files:

ins_tr_20.tsv, ins_not_tr_20.tsv, del_tr_20.tsv, del_not_tr_20.tsv,
ins_tr_50.tsv, ins_not_tr_50.tsv, del_tr_50.tsv, del_not_tr_50.tsv

Every file is a TSV with format `SVTYPE \t SVLEN \t NumNeighbors`. A first line
that is not numeric in the last column is treated as a header and skipped.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Column of the TSV
COL_NUM_NEIGHBORS = 2

# Colors: (INS, DEL)
COLORS_NOT_TR = ("#b01010", "#1050b0")
COLORS_TR = ("#f4a0a0", "#9cc4ee")

FIXED_BINS = [1, 2, 4, 8, 16, 20]
INCREMENT = 10
MAX_N = 10000


def load(path):
    values = []
    with open(path) as f:
        for i, line in enumerate(f):
            tokens = line.rstrip("\n").split("\t")
            if len(tokens) <= COL_NUM_NEIGHBORS:
                continue
            try:
                values.append(int(tokens[COL_NUM_NEIGHBORS]))
            except ValueError:
                if i == 0:
                    continue  # Header
                raise
    return np.array(values, dtype=np.int64)


def get_bins():
    bins = list(FIXED_BINS)
    while bins[-1] < MAX_N:
        bins.append(bins[-1] + INCREMENT)
    return np.array(bins, dtype=np.int64)


def plot_curve(ax, values, bins, label, color, linestyle):
    n = max(len(values), 1)
    fractions = [np.count_nonzero(values >= b) / n for b in bins]
    ax.plot(bins, fractions, color=color, linestyle=linestyle, linewidth=2, label=label)


def main():
    output_png = sys.argv[1]
    title = sys.argv[2]
    indir = sys.argv[3] if len(sys.argv) > 3 else "."

    # (file name, label, color, linestyle)
    curves = []
    for length, linestyle in (("20", "--"), ("50", "-")):
        for svtype, idx in (("ins", 0), ("del", 1)):
            for tr, colors in (("not_tr", COLORS_NOT_TR), ("tr", COLORS_TR)):
                name = f"{svtype}_{tr}_{length}.tsv"
                label = rf"$\geq${length}bp, {'non-TR' if tr == 'not_tr' else 'TR'}, {svtype.upper()}"
                curves.append((name, label, colors[idx], linestyle))
    data = [load(os.path.join(indir, c[0])) for c in curves]
    bins = get_bins()

    fig, ax = plt.subplots(figsize=(10, 8))
    for (name, label, color, linestyle), values in zip(curves, data):
        plot_curve(ax, values, bins, label, color, linestyle)
    ax.set_title(title)
    ax.set_xscale("log")
    ax.set_xlim(1, bins[-1])
    ax.set_ylim(0, 1)
    ax.set_box_aspect(1)  # Square axes
    ax.set_axisbelow(True)
    ax.grid(True, color="lightgray", linewidth=0.8)
    ax.set_xlabel("Number of neighbors")
    ax.set_ylabel(r"Fraction of calls with $\geq$ n neighbors")
    ax.legend(loc="upper right", fontsize="small")
    fig.tight_layout()
    fig.savefig(output_png, dpi=150)


if __name__ == "__main__":
    main()
