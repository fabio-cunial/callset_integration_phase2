#!/usr/bin/env python3
"""
Plots histograms of the smallest distances computed by `SmallestDistance.java`:
- INS to a call of any type;
- INS to INS;
- DEL to a call of any type;
- DEL to DEL.

Calls with no call of the given type within the max distance (value -1) are
shown as a separate bar to the left of zero, and their number is reported in
each title.

Usage: plot_smallest_distance.py output.png n_bins input1.bed [input2.bed ...]

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


def plot(ax, distances, n_bins, title):
    found = distances[distances >= 0]
    n_missing = len(distances) - len(found)
    if len(found) > 0:
        edges = np.histogram_bin_edges(found, bins=n_bins)
        ax.hist(found, bins=edges, color="steelblue", label="found")
        width = edges[1] - edges[0] if edges[1] > edges[0] else 1
    else:
        edges, width = np.array([0]), 1
    # Calls with distance -1 are drawn as a separate bar to the left of zero.
    ax.bar(edges[0] - width, n_missing, width=width, align="edge", color="darkorange", label="not found within 1Mbp")
    ax.legend()
    ax.set_axisbelow(True)
    ax.grid(True, color="lightgray", linewidth=0.8)
    ax.set_box_aspect(1)
    ax.set_title(f"{title}\n(found={len(found)}, not found={n_missing})")
    ax.set_xlabel("Smallest distance (bp)")
    ax.set_ylabel("Number of calls")


def main():
    output_png = sys.argv[1]
    n_bins = int(sys.argv[2])
    data = load(sys.argv[3:])
    svtype, to_any, to_del, to_ins = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
    is_ins = svtype == TYPE_INS
    is_del = svtype == TYPE_DEL

    fig, axes = plt.subplots(2, 2, figsize=(11, 11))
    plot(axes[0, 0], to_any[is_ins], n_bins, "INS -> any")
    plot(axes[0, 1], to_ins[is_ins], n_bins, "INS -> INS")
    plot(axes[1, 0], to_any[is_del], n_bins, "DEL -> any")
    plot(axes[1, 1], to_del[is_del], n_bins, "DEL -> DEL")
    fig.tight_layout()
    fig.savefig(output_png, dpi=150)


if __name__ == "__main__":
    main()
