#!/usr/bin/env python3
"""
Plots two 2D-histogram heatmaps comparing old vs. new values.

Input CSV columns (no header assumed by default):
    ID,CAL_SENS_OLD,SCORE_OLD,CAL_SENS_NEW,SCORE_NEW

Left panel : heatmap of CAL_SENS_OLD (x) vs. CAL_SENS_NEW (y).
Right panel: heatmap of SCORE_OLD  (x) vs. SCORE_NEW  (y).
Each cell holds the number of records falling into that (x,y) bin.

Usage:
    python3 plot_old_new_heatmaps.py <input.csv> [n_bins] [has_header]

    n_bins     : number of bins per axis (default 100)
    has_header : 1 if the first CSV line is a header, else 0 (default 1)
"""

import sys
import csv

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(__doc__)
        sys.exit(1)

    input_csv = sys.argv[1]
    n_bins = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    has_header = (int(sys.argv[3]) if len(sys.argv) > 3 else 1) != 0

    cal_old, cal_new, score_old, score_new = [], [], [], []
    with open(input_csv, newline="") as f:
        reader = csv.reader(f)
        if has_header:
            next(reader, None)
        for row in reader:
            if len(row) < 5:
                continue
            try:
                co = float(row[1])
                so = float(row[2])
                cn = float(row[3])
                sn = float(row[4])
            except ValueError:
                # Skip malformed / non-numeric lines.
                continue
            cal_old.append(co)
            score_old.append(so)
            cal_new.append(cn)
            score_new.append(sn)

    cal_old = np.asarray(cal_old)
    cal_new = np.asarray(cal_new)
    score_old = np.asarray(score_old)
    score_new = np.asarray(score_new)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    panels = [
        (axes[0], cal_old, cal_new, "CAL_SENS_OLD", "CAL_SENS_NEW", ""),
        (axes[1], score_old, score_new, "SCORE_OLD", "SCORE_NEW", ""),
    ]

    for ax, x, y, xlabel, ylabel, title in panels:
        r2 = np.corrcoef(x, y)[0, 1] ** 2 if x.size > 1 else float("nan")
        title = "{} ($R^2$ = {:.3f})".format(title, r2).strip()
        counts, xedges, yedges = np.histogram2d(x, y, bins=n_bins)
        # histogram2d indexes [x, y]; transpose so rows=y, cols=x for imshow.
        im = ax.imshow(
            counts.T,
            origin="lower",
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            aspect="auto",
            cmap="viridis",
            norm=LogNorm(vmin=1),
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, color="lightgray", alpha=1, linewidth=1)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Number of records")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
