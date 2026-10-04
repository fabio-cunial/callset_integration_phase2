#!/usr/bin/env python3
"""
Python port of SV_Integration_RegenotypingAnalysis_plot.m

Draws one figure per evaluation threshold (50bp, 20bp). Each figure has:
  row 1 = precision / recall / GT concordance
  row 2 = Mendelian error rate
  row 3 = de novo rate
  column 1 = inside TRs, column 2 = outside TRs

Usage:  python3 SV_Integration_RegenotypingAnalysis_plot.py [input_dir]
The CSVs are read from input_dir (default: the directory of this script).
"""

import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MIN_N_SAMPLES = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
LABELS = ["T", "1", "2", "4", "8", "16", "32", "64", "128", "256", "512", "1024", "2048"]
DELTA = 0.4
EVAL_THRESHOLDS = ["50bp", "20bp"]
REGIONS = [("tr", "TR regions"), ("not_tr", "non-TR regions")]
# Call-length range and chromosome covered by each threshold's CSVs.
TITLE_PREFIXES = {"20bp": "[20..50) bp, chr6", "50bp": "[50..10001) bp, chr6"}

# Okabe-Ito categorical hues: distinguishable under protan/deutan/tritan vision,
# unlike the matplotlib defaults (whose green and orange collapse for protanopes).
BLUE = "#0072B2"
VERMILLION = "#D55E00"
GREEN = "#009E73"
# Markers overlap heavily at every x, so draw them small and semi-transparent:
# overlap then reads as density instead of as one opaque blob.
# The circles stay hollow on purpose: that is what separates them from the dots
# once both wear the same hue, and overlapping rings stay readable where discs
# would merge into a blob.
DOT = dict(linestyle="none", marker=".", markersize=6.5, alpha=0.8)
CIRCLE = dict(linestyle="none", marker="o", markersize=8, markerfacecolor="none",
              markeredgewidth=1.1, alpha=0.8)

# Deterministic jitter, so repeated runs give identical figures.
RNG = np.random.default_rng(0)

INPUT_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))


def load(name):
    """Loads a comma-separated numeric matrix, tolerating a trailing comma."""
    rows = []
    with open(os.path.join(INPUT_DIR, name)) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append([float(x) for x in line.split(",") if x != ""])
    return np.array(rows, dtype=float)


def jitter(center, n):
    """n x-coordinates scattered uniformly in [center-DELTA/2, center+DELTA/2)."""
    return center - DELTA / 2 + RNG.random(n) * DELTA


def ratio(numerator, denominator):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(numerator, denominator,
                         out=np.full_like(numerator, np.nan, dtype=float),
                         where=denominator != 0)


def decorate(ax, title, ylabel, ylim):
    ax.set_xticks(range(1, len(LABELS) + 1))
    ax.set_xticklabels(LABELS, rotation=45)
    ax.set_xlabel("Min n. samples")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, color="0.85", linewidth=0.8)
    ax.set_xlim(0, len(LABELS) + 1)
    ax.set_ylim(*ylim)
    ax.set_box_aspect(1)


def dot(color):
    return Line2D([], [], color=color, markeredgecolor=color,
                  **dict(DOT, alpha=1, markersize=7))


def circle(color):
    return Line2D([], [], color=color, markeredgecolor=color,
                  **dict(CIRCLE, alpha=1, markersize=7))


# ----------------------------- Precision/recall -------------------------------

def plot_precision_recall(ax, threshold, region, title):
    # Column 0 = precision, 1 = recall, 3 = GT concordance.
    measures = ((0, BLUE), (1, VERMILLION), (3, GREEN))
    for coverage, style in (("15x", DOT), ("30x", CIRCLE)):
        A = load("precision_recall_%s_%s_%s.csv" % (coverage, threshold, region))
        for i in range(A.shape[0]):
            # Truvari collapse.
            X = jitter(1, 1)
            for column, color in measures:
                ax.plot(X, [A[i, column]], color=color, markeredgecolor=color, **style)
            # Re-genotyping.
            n = len(MIN_N_SAMPLES)
            X = 1 + jitter(0, n) + np.arange(1, n + 1)
            for column, color in measures:
                Y = [A[i, j * 4 + column] for j in range(1, n + 1)]
                ax.plot(X, Y, color=color, markeredgecolor=color, **style)
        if coverage == "15x":
            # Baseline: the mean over the 15x truvari-collapse points at "T", so
            # the later x values can be read against the starting point.
            for column, color in measures:
                ax.axhline(A[:, column].mean(), xmin=1 / (len(LABELS) + 1),
                           xmax=len(LABELS) / (len(LABELS) + 1),
                           linestyle="--", color=color, linewidth=1, zorder=1)

    decorate(ax, title, " ", (0.4, 1))
    ax.legend([dot(BLUE), dot(VERMILLION), dot(GREEN),
               circle(BLUE), circle(VERMILLION), circle(GREEN)],
              ["15x precision", "15x recall", "15x GT conc.",
               "30x precision", "30x recall", "30x GT conc."],
              loc="lower left", ncol=2, framealpha=0.9)


# --------------------- Mendelian error / de novo rate -------------------------

def plot_rate(ax, prefix, threshold, region, title, ylabel, step, rate_at):
    """Mendelian error and de novo rate share a layout; only `step` (the number
    of CSV columns per x position) and `rate_at` (the rate formula) differ."""
    series = []
    for coverage, cohort, color in (("15x", "control", BLUE),
                                    ("15x", "aou", BLUE),
                                    ("30x", "aou", VERMILLION)):
        for suffix, style, missing_to_ref in (("", DOT, False),
                                              ("_no_missing", CIRCLE, True)):
            M = load("%s_%s_%s_%s_%s%s.csv"
                     % (prefix, coverage, threshold, cohort, region, suffix))
            series.append((M, style, color, missing_to_ref))

    ncolumns = series[2][0].shape[1]  # The 15x AoU matrix sets the loop bounds.
    for i in range(step, ncolumns + 1, step):  # 1-based index, as in the .m script.
        for M, style, color, _ in series:
            if M.shape[1] < i:
                continue
            ax.plot(jitter(i / step, M.shape[0]), rate_at(M, i),
                    color=color, markeredgecolor=color, **style)

    # Baselines: the mean of the "missing->ref" points at "T", per coverage, so
    # the later x values can be read against the starting point.
    for color in (BLUE, VERMILLION):
        Y = np.concatenate([rate_at(M, step) for M, _, c, missing_to_ref in series
                            if missing_to_ref and c == color])
        ax.axhline(np.nanmean(Y), xmin=1 / (len(LABELS) + 1),
                   xmax=len(LABELS) / (len(LABELS) + 1),
                   linestyle="--", color=color, linewidth=1, zorder=1)

    decorate(ax, title, ylabel, (0, 0.25))
    ax.legend([dot(BLUE), circle(BLUE), dot(VERMILLION), circle(VERMILLION)],
              ["15x (control and AoU)", "15x, missing->ref",
               "30x AoU", "30x AoU, missing->ref"],
              loc="upper right", ncol=1, framealpha=0.9)


def plot_mendelian_error(ax, threshold, region, title):
    plot_rate(ax, "mendelian_error", threshold, region, title, "Mendelian error rate",
              2, lambda M, i: ratio(M[:, i - 1], M[:, i - 1] + M[:, i - 2]))


def plot_denovo(ax, threshold, region, title):
    plot_rate(ax, "denovo", threshold, region, title, "De novo rate",
              6, lambda M, i: ratio(M[:, i - 2], M[:, i - 1]))


def build_figure(threshold):
    # 5x5 inches per subplot, 3 rows x 2 columns.
    fig, axes = plt.subplots(3, 2, figsize=(5 * len(REGIONS), 5 * 3))
    for column, (region, region_title) in enumerate(REGIONS):
        title = "%s, %s" % (TITLE_PREFIXES[threshold], region_title)
        plot_precision_recall(axes[0][column], threshold, region, title)
        plot_mendelian_error(axes[1][column], threshold, region, title)
        plot_denovo(axes[2][column], threshold, region, title)
    # The x axis is the same in every row: label it only on the bottom one.
    for row in axes[:-1]:
        for ax in row:
            ax.set_xlabel("")
    # The y axis is the same in every column: label it only on the leftmost one.
    for row in axes:
        for ax in row[1:]:
            ax.set_ylabel("")
    fig.tight_layout()
    return fig


def main():
    for threshold in EVAL_THRESHOLDS:
        fig = build_figure(threshold)
        out = os.path.join(INPUT_DIR, "regenotyping_analysis_%s.png" % threshold)
        fig.savefig(out, dpi=150)
        print("Wrote %s" % out)
    if matplotlib.get_backend().lower() != "agg":
        plt.show()


if __name__ == "__main__":
    main()
