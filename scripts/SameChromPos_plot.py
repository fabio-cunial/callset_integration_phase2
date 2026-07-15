#!/usr/bin/env python3
"""
Plots ratios of the fields produced by SameChromPos.java.

Input: a CSV file where every row has the following fields:
0: total number of records;
1: number of records with the same CHROM,POS as another record;

2: number of distinct CHROM,POS pairs;
3: number of CHROM,POS pairs with an INS and a DEL;
4: number of CHROM,POS pairs with multiple INS;
5: number of CHROM,POS pairs with multiple DEL;
6: number of CHROM,POS pairs with >1 record.

A second input text file, with one integer per row, is shown as a histogram
in a second panel of the same figure.

A third input CSV file, with two floats per row, is shown as two histograms
in a third and fourth panel of the same figure.

A fourth input CSV file, with three floats per row, contributes an extra
coordinate to the first panel: for every row, the third column divided by the
second column is plotted as a jittered circle with a box-and-whiskers overlay,
just like the other coordinates.

Usage:
  ./plot_same_chrom_pos.py <input.csv> <cluster_sizes.txt> <deltas.csv> <trs.csv> [title]
"""

import sys
import csv
import random
import numpy as np
import matplotlib.pyplot as plt


RATIOS = [
    ("Records with \n same (CHR,POS) \n as another record \n (fraction of all records)", 1, 0),
    ("Problematic (CHR,POS) \n pairs (fraction of all \n distinct pairs)", 6, 2),
    ("(CHR,POS) pairs \n with INS and DEL \n (fraction of \n all problematic pairs)", 3, 6),
    ("(CHR,POS) pairs with \n >=2 INS (fraction \n of all problematic pairs)", 4, 6),
    ("(CHR,POS) pairs with \n >=2 DEL (fraction \n of all problematic \n pairs)", 5, 6),
]


def main():
    if len(sys.argv) < 5:
        sys.stderr.write("Usage: %s <input.csv> <cluster_sizes.txt> <deltas.csv> <trs.csv> [title]\n" % sys.argv[0])
        sys.exit(1)
    input_csv = sys.argv[1]
    cluster_sizes_txt = sys.argv[2]
    deltas_csv = sys.argv[3]
    trs_csv = sys.argv[4]
    title = sys.argv[5] if len(sys.argv) > 5 else " "

    rows = []
    with open(input_csv, newline="") as f:
        for record in csv.reader(f):
            record = [c.strip() for c in record if c.strip() != ""]
            if not record:
                continue
            fields = [float(c) for c in record[:7]]
            if len(fields) < 7:
                sys.stderr.write("Skipping row with fewer than 7 fields: %s\n" % record)
                continue
            rows.append(fields)

    if not rows:
        sys.stderr.write("No valid rows found in %s\n" % input_csv)
        sys.exit(1)

    # Read the integers for the second panel, one per row.
    integers = []
    with open(cluster_sizes_txt) as f:
        for line in f:
            line = line.strip()
            if line == "":
                continue
            integers.append(int(line))

    if not integers:
        sys.stderr.write("No valid integers found in %s\n" % cluster_sizes_txt)
        sys.exit(1)

    # Read the (float, float) pairs for the third and fourth panels.
    first_values = []
    second_values = []
    with open(deltas_csv, newline="") as f:
        for record in csv.reader(f):
            record = [c.strip() for c in record if c.strip() != ""]
            if not record:
                continue
            if len(record) < 2:
                sys.stderr.write("Skipping row with fewer than 2 fields: %s\n" % record)
                continue
            first_values.append(float(record[0]))
            second_values.append(float(record[1]))

    if not first_values:
        sys.stderr.write("No valid rows found in %s\n" % deltas_csv)
        sys.exit(1)

    # Read the (float, float, float) triples for the extra coordinate of the
    # first panel: the ratio is the third column divided by the second column.
    extra_ratios = []
    with open(trs_csv, newline="") as f:
        for record in csv.reader(f):
            record = [c.strip() for c in record if c.strip() != ""]
            if not record:
                continue
            if len(record) < 3:
                sys.stderr.write("Skipping row with fewer than 3 fields: %s\n" % record)
                continue
            second = float(record[1])
            if second == 0:
                continue
            extra_ratios.append(float(record[2]) / second)

    if not extra_ratios:
        sys.stderr.write("No valid rows found in %s\n" % trs_csv)
        sys.exit(1)

    # Deterministic jitter for reproducible plots.
    random.seed(0)
    jitter_width = 0.15

    # Label for the extra coordinate contributed by <trs.csv>.
    extra_label = "(CHR,POS) pairs in TRs \n (fraction of all \n problematic pairs)"

    # Collect the finite ratio values per coordinate for the box-and-whiskers.
    # The extra coordinate from <extra.csv> is appended after the RATIOS ones.
    values_per_x = []
    for x, (label, num, den) in enumerate(RATIOS):
        values = []
        for fields in rows:
            if fields[den] == 0:
                continue
            values.append(fields[num] / fields[den])
        values_per_x.append(values)
    values_per_x.append(list(extra_ratios))

    n_x = len(RATIOS) + 1
    x_labels = [label for label, _, _ in RATIOS] + [extra_label]

    fig, ((ax, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 12))

    # Box-and-whiskers underlay, one box per coordinate.
    ax.boxplot(values_per_x, positions=range(n_x), widths=0.5,
               showfliers=False, zorder=1,
               medianprops=dict(color="tab:red"),
               boxprops=dict(color="gray"),
               whiskerprops=dict(color="gray"),
               capprops=dict(color="gray"))

    # Jittered circles, one per row, overlaid on top of the boxes.
    for x, (label, num, den) in enumerate(RATIOS):
        for fields in rows:
            ratio = fields[num] / fields[den] if fields[den] != 0 else float("nan")
            xj = x + random.uniform(-jitter_width, jitter_width)
            ax.plot(xj, ratio, marker="o", markersize=8, markerfacecolor="none",
                    markeredgecolor="tab:blue", linestyle="none", alpha=0.6, zorder=2)

    # Jittered circles for the extra coordinate from <extra.csv>.
    for ratio in extra_ratios:
        xj = len(RATIOS) + random.uniform(-jitter_width, jitter_width)
        ax.plot(xj, ratio, marker="o", markersize=8, markerfacecolor="none",
                markeredgecolor="tab:blue", linestyle="none", alpha=0.6, zorder=2)

    ax.set_xticks(range(n_x))
    ax.set_xticklabels(x_labels)
    ax.set_xlim(-0.5, n_x - 0.5)
    #ax.set_ylim(top=0.07)
    ax.set_ylabel("Fraction")
    ax.set_title("%s (%d HPRC+HGSVC samples)" % (title, len(rows)))
    ax.grid(True, axis="y", linestyle=":", alpha=0.5)

    # Second panel: histogram of cluster sizes.
    ax2.hist(integers, bins="auto", color="tab:blue", edgecolor="black", alpha=0.7)
    ax2.set_xlabel("Number of records of a (CHR,POS) pair")
    ax2.set_ylabel("Number of problematic (CHR,POS) pairs")
    ax2.set_title("%s, number of records of problematic pairs, %d (CHROM,POS) pairs" % (title, len(integers)))
    ax2.grid(True, axis="y", linestyle=":", alpha=0.5)

    # Weights that make each histogram's bar heights sum to 1 (fraction of mass).
    first_weights = [1.0 / len(first_values)] * len(first_values)
    second_weights = [1.0 / len(second_values)] * len(second_values)

    # Auto bin edges (unweighted), so the weighted histograms keep automatic binning.
    first_bins = np.histogram_bin_edges(first_values, bins="auto")
    second_bins = np.histogram_bin_edges(second_values, bins="auto")

    # Third panel: histogram of maxLen-minLen.
    ax3.hist(first_values, bins=first_bins, weights=first_weights,
             color="tab:green", edgecolor="black", alpha=0.7)
    ax3.set_xlabel("maxLen-minLen")
    ax3.set_ylabel("Fraction of problematic (CHR,POS) pairs")
    ax3.set_title("%s, maxLen-minLen, %d (CHR,POS) pairs" % (title, len(first_values)))
    ax3.grid(True, axis="y", linestyle=":", alpha=0.5)

    # Fourth panel: histogram of the relative deltas.
    ax4.hist(second_values, bins=second_bins, weights=second_weights,
             color="tab:orange", edgecolor="black", alpha=0.7)
    ax4.set_xlabel("(maxLen-minLen)/minLen")
    ax4.set_ylabel("Fraction of problematic (CHR,POS) pairs")
    ax4.set_title("%s, (maxLen-minLen)/minLen, %d (CHR,POS) pairs" % (title, len(second_values)))
    ax4.grid(True, axis="y", linestyle=":", alpha=0.5)

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
