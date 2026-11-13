#!/usr/bin/env python3
"""
Draws a figure with four heatmaps. Rows: top=SVTYPE=DEL, bottom=SVTYPE=INS.
Columns: left=calls inside tandem repeats, right=calls outside tandem repeats.
In each heatmap the X axis is SVLEN, the Y axis is ALLELE_FREQUENCY, and cell
(x,y) contains the number of rows with that combination of values. Every
heatmap of the top row has, glued to its right, a narrow histogram of the
ALLELE_FREQUENCY values of the same calls (the fraction of the panel's calls in
each AF bin, scaled to its own largest bar).

Input: two headerless TSVs (inside TRs, outside TRs), each with exactly the
columns SVTYPE, SVLEN, AC, in this order. SVLEN is used in absolute value, so
negative lengths for deletions are fine. ALLELE_FREQUENCY is computed as
AC/(2*n_samples), i.e. missing genotypes count as reference. The Y axis is
always log10, so records with AC=0 are discarded.

Binning happens in integer space on both axes: the AC and SVLEN edges are
rounded to integers and de-duplicated, and only the AC edges are then rescaled
by 2*n_samples for drawing. Every bin therefore covers a whole number of AC
(resp. bp) values, which removes the empty/doubled stripes you otherwise get at
low AF (resp. short SVLEN) where a log bin would be narrower than one allele
(resp. one bp). Hence x_bins and y_bins are targets: the actual number of bins
can be smaller.

Usage: AfLen.py <inside_tr.tsv[.gz]> <outside_tr.tsv[.gz]> <output.png> <n_samples> [x_bins] [y_bins] [log_x]
  n_samples: total number of samples in the cohort
  x_bins: number of SVLEN bins (default 100)
  y_bins: number of ALLELE_FREQUENCY bins (default 50)
  log_x:  1 (default) for a log10 X axis, 0 for linear
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

TYPES = ["DEL", "INS"]  # top to bottom
REGIONS = ["TR regions", "non-TR regions"]  # left to right

# Figure geometry, in inches.
PANEL = 4.0          # side of a square heatmap
MARGINAL = 0.5       # width of a marginal histogram
MARGINAL_GAP = 0.04  # between a heatmap and its marginal
GUTTER = 0.9         # between the two columns; holds the right y tick labels/label
ROW_GAP = 0.9        # between the two rows; holds the bottom-row titles
MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM, MARGIN_TOP = 0.9, 0.1, 0.7, 0.4

MARGINAL_COLOR = "#2a78d6"


def load(path, log_x):
    df = pd.read_csv(path, sep="\t", header=None, names=["SVTYPE", "SVLEN", "AC"])
    df["SVLEN"] = pd.to_numeric(df["SVLEN"], errors="coerce").abs()
    df["AC"] = pd.to_numeric(df["AC"], errors="coerce")
    df = df.dropna()
    df = df[df["AC"] >= 1]
    df["AC"] = df["AC"].astype(np.int64)
    df["SVLEN"] = df["SVLEN"].astype(np.int64)
    if log_x:
        df = df[df["SVLEN"] > 0]
    return df[df["SVTYPE"].isin(TYPES)]


def log_integer_edges(lo, hi, n_bins):
    """Log-spaced bin edges over the integers [lo,hi], snapped to integers and
    de-duplicated so no bin is narrower than one unit.

    Returns half-open edges e with e[0]=lo and e[-1]=hi+1, i.e. bin k holds the
    values in [e[k], e[k+1]-1].
    """
    raw = np.geomspace(lo, hi + 1, n_bins + 1)
    edges = np.unique(np.round(raw).astype(np.int64))
    edges[0], edges[-1] = lo, hi + 1
    return np.unique(edges)


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    in_tsv, out_tsv, out_png = sys.argv[1], sys.argv[2], sys.argv[3]
    n_samples = int(sys.argv[4])
    if n_samples <= 0:
        sys.exit("n_samples must be positive.")
    x_bins = int(sys.argv[5]) if len(sys.argv) > 5 else 100
    y_bins = int(sys.argv[6]) if len(sys.argv) > 6 else 50
    log_x = (sys.argv[7] != "0") if len(sys.argv) > 7 else True

    n_alleles = 2 * n_samples  # AF = AC / (2N)

    dfs = [load(in_tsv, log_x), load(out_tsv, log_x)]
    all_df = pd.concat(dfs)
    if all_df.empty:
        sys.exit("No INS/DEL rows found.")
    # A wrong N (or a non-diploid callset) shows up here rather than as a
    # silently mislabelled axis.
    if all_df["AC"].max() > n_alleles:
        sys.exit(f"AC up to {all_df['AC'].max():,} exceeds 2N = {n_alleles:,} for N = {n_samples:,} samples.")

    # Bin edges are shared by all panels so that they are comparable.
    # SVLEN is binned in integer space too, so that short calls, where a log
    # bin would be narrower than 1bp, do not produce empty/doubled stripes.
    lo, hi = int(all_df["SVLEN"].min()), int(all_df["SVLEN"].max())
    if log_x:
        len_edges = log_integer_edges(lo, hi, x_bins)
    else:
        len_edges = np.unique(np.round(np.linspace(lo, hi + 1, x_bins + 1)).astype(np.int64))
    # Integer length k occupies [k-0.5, k+0.5) on the drawn axis.
    x_edges = len_edges - 0.5
    # AC bins span every possible nonzero AC, i.e. AF from 1/(2N) to 1. They are
    # rescaled to AF only for drawing, so integer AC k occupies
    # [(k-0.5)/2N, (k+0.5)/2N).
    ac_edges = log_integer_edges(1, n_alleles, y_bins)
    y_edges = (ac_edges - 0.5) / n_alleles

    hists = {}
    for r, df in enumerate(dfs):
        for t_idx, t in enumerate(TYPES):
            sub = df[df["SVTYPE"] == t]
            h, _, _ = np.histogram2d(sub["SVLEN"], sub["AC"], bins=[len_edges, ac_edges])
            hists[(t_idx, r)] = h.T  # rows = AF, columns = SVLEN

    vmax = max(h.max() for h in hists.values())
    norm = LogNorm(vmin=1, vmax=max(vmax, 2))  # log color scale; empty cells are masked

    # Axes are laid out by hand, in inches, so that every heatmap stays exactly
    # square while its marginal stays glued to it at a fixed width. Both columns
    # reserve room for a marginal, so the bottom heatmaps align with the top ones.
    column_width = PANEL + MARGINAL_GAP + MARGINAL
    col_x = (MARGIN_LEFT, MARGIN_LEFT + column_width + GUTTER)
    row_y = (MARGIN_BOTTOM + PANEL + ROW_GAP, MARGIN_BOTTOM)  # top row, bottom row
    fig_w = col_x[1] + column_width + MARGIN_RIGHT
    fig_h = MARGIN_BOTTOM + 2 * PANEL + ROW_GAP + MARGIN_TOP
    fig = plt.figure(figsize=(fig_w, fig_h))

    def add_axes(x, y, w, h, **kwargs):
        return fig.add_axes((x / fig_w, y / fig_h, w / fig_w, h / fig_h), **kwargs)

    first = None
    for (t_idx, r), h in hists.items():
        ax = add_axes(col_x[r], row_y[t_idx], PANEL, PANEL, sharex=first, sharey=first)
        first = first or ax
        ax.pcolormesh(x_edges, y_edges, np.ma.masked_equal(h, 0), norm=norm, cmap="viridis")
        if log_x:
            ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{TYPES[t_idx]}, {REGIONS[r]}.") #ax.set_title(f"{TYPES[t_idx]}, {REGIONS[r]} (n={int(h.sum()):,})")
        if t_idx == 1:
            ax.set_xlabel("Length (bp)")
        if r == 0:
            ax.set_ylabel("AF")
        if t_idx == 0:
            # Marginal: the fraction of the panel's calls in each AF bin. Summing
            # the heatmap's columns rather than re-histogramming keeps it the exact
            # row marginal of the heatmap. It has no axis and is scaled to its own
            # largest bar, since it is read only for its shape. Bar lengths are
            # log10: every bar starts at the fraction of half a call, so a bin with
            # one call is still visible and an empty bin has no bar.
            totals = h.sum(axis=1)
            if totals.sum():
                fractions = totals / totals.sum()
                floor = 0.5 / totals.sum()
                ax_dist = add_axes(col_x[r] + PANEL + MARGINAL_GAP, row_y[t_idx], MARGINAL, PANEL, sharey=ax)
                ax_dist.barh(y_edges[:-1], np.where(fractions > 0, fractions - floor, 0), left=floor,
                             height=np.diff(y_edges), align="edge",
                             color=MARGINAL_COLOR, edgecolor="white", linewidth=0.8)
                ax_dist.set_xscale("log")
                ax_dist.set_xlim(floor, fractions.max() * 1.05)
                ax_dist.axis("off")
    fig.savefig(out_png, dpi=200, bbox_inches="tight", pad_inches=0.01)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
