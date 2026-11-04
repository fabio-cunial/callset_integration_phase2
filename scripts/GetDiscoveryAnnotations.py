#!/usr/bin/env python3
"""Plot the SUPP and SCORE histograms of the first records of a TSV produced by
`GetDiscoveryAnnotations.java` (`ID \t SUPP_HISTOGRAM \t SCORE_HISTOGRAM`).

Every record is one row of two side-by-side heatmaps that share the y axis:
- left: SUPP histogram, one column per nonzero bitmap value (PBSV=1,
  SNIFFLES=2, PAV=4). The no-support bin (0) is not plotted;
- right: SCORE histogram, one column per bin in [SCORE_MIN..SCORE_MAX]. The
  first and last bins (scores outside the range) are not plotted.

Every cell shows the raw number of samples in that bin. Color is log-scaled
and shared by both plots; empty cells are left blank.

Usage:
    PlotDiscoveryAnnotations.py INPUT_TSV SCORE_MIN SCORE_MAX OUTPUT_PNG [N_LINES=100]
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

SUPP_LABELS = ["pbsv", "snf", "pbsv+snf", "pav", "pbsv+pav", "snf+pav", "all"]
SEQUENTIAL = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
SURFACE = "#ffffff"
INK_SECONDARY = "#55554f"


def load(path, n_lines):
    ids, supp, score = [], [], []
    with open(path) as f:
        for line in f:
            if len(ids) == n_lines:
                break
            id, s, c = line.rstrip("\n").split("\t")
            ids.append(id)
            supp.append([int(x) for x in s.split(",")])
            score.append([int(x) for x in c.split(",")])
    return ids, np.array(supp, dtype=float), np.array(score, dtype=float)


def main():
    input_tsv = sys.argv[1]
    score_min = float(sys.argv[2])
    score_max = float(sys.argv[3])
    output_png = sys.argv[4]
    n_lines = int(sys.argv[5]) if len(sys.argv) > 5 else 100

    ids, supp, score = load(input_tsv, n_lines)
    supp = supp[:, 1:]  # Dropping the no-support bin
    score = score[:, 1:-1]  # Dropping out-of-range bins
    n_rows, n_bins = len(ids), score.shape[1]
    supp = np.ma.masked_equal(supp, 0)
    score = np.ma.masked_equal(score, 0)

    norm = LogNorm(vmin=1, vmax=max(supp.max(), score.max(), 2))
    cmap = LinearSegmentedColormap.from_list("seq", SEQUENTIAL)
    cmap.set_bad(SURFACE)

    height = max(4.0, 0.12 * n_rows + 1.5)
    fig, (ax_supp, ax_score) = plt.subplots(
        1, 2, figsize=(12, height), sharey=True,
        gridspec_kw={"width_ratios": [len(SUPP_LABELS), max(n_bins, len(SUPP_LABELS))], "wspace": 0.04},
    )
    for ax, counts in ((ax_supp, supp), (ax_score, score)):
        image = ax.imshow(counts, aspect="auto", cmap=cmap, norm=norm,
                          interpolation="nearest")
        # 1px surface gaps between cells, so adjacent rows are separable.
        ax.set_xticks(np.arange(-0.5, counts.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
        ax.grid(which="minor", color=SURFACE, linewidth=0.6)
        ax.tick_params(which="minor", length=0)
        ax.tick_params(colors=INK_SECONDARY, labelsize=7)
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        for spine in ax.spines.values():
            spine.set_visible(False)

    ax_supp.set_xticks(range(len(SUPP_LABELS)))
    ax_supp.set_xticklabels(SUPP_LABELS, rotation=60, ha="left", rotation_mode="anchor")
    ax_supp.set_title("SUPP", fontsize=10, pad=45)

    edges = np.linspace(score_min, score_max, n_bins + 1)
    step = max(1, n_bins // 10)
    ticks = np.arange(0, n_bins + 1, step)
    ax_score.set_xticks(ticks - 0.5)
    ax_score.set_xticklabels([f"{edges[i]:g}" for i in ticks])
    ax_score.set_title(f"SCORE in [{score_min:g}..{score_max:g}] (out-of-range bins omitted)",
                       fontsize=10, pad=45)
    ax_score.set_xlabel("score", fontsize=8, color=INK_SECONDARY)
    ax_score.tick_params(axis="y", which="major", left=False)

    if n_rows <= 100:
        ax_supp.set_yticks(range(n_rows))
        ax_supp.set_yticklabels(ids, fontsize=max(3, min(7, 400 // n_rows)))
    ax_supp.set_ylabel("record", fontsize=8, color=INK_SECONDARY)

    cbar = fig.colorbar(image, ax=[ax_supp, ax_score], fraction=0.02, pad=0.01)
    cbar.set_label("number of samples (log)", fontsize=8, color=INK_SECONDARY)
    cbar.ax.tick_params(labelsize=7, colors=INK_SECONDARY)
    cbar.outline.set_visible(False)

    fig.savefig(output_png, dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
