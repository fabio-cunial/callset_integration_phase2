#!/usr/bin/env python3
"""
Plots one histogram per INFO field of the TSV produced by
wdl/InvestigateComplexDuplicationsPrime.wdl, all as panels of a single figure.

Columns (no header):
ID, SVLEN, TruScore, PctSeqSimilarity, PctSizeSimilarity, PctRecOverlap,
StartDistance, EndDistance, SizeDiff

Usage: InvestigateComplexDuplicationsPrime.py <input.tsv> <output.png> [n_bins]
"""
import sys
import math

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLUMNS = ["ID", "SVLEN", "TruScore", "PctSeqSimilarity", "PctSizeSimilarity",
           "PctRecOverlap", "StartDistance", "EndDistance", "SizeDiff"]


LOG_FIELDS = {"SVLEN", "StartDistance", "SizeDiff"}


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    in_tsv, out_png = sys.argv[1], sys.argv[2]
    n_bins = int(sys.argv[3]) if len(sys.argv) > 3 else 50

    # bcftools prints missing values as '.'
    df = pd.read_csv(in_tsv, sep="\t", header=None, names=COLUMNS,
                     na_values=["."])
# EndDistance is skipped, since it seems identical to StartDistance;
    # TruScore, since it seems identical to PctSeqSimilarity.
    fields = [c for c in COLUMNS[1:] if c not in ("EndDistance", "TruScore")]

    n_cols = 3
    n_rows = math.ceil((len(fields) + 1) / n_cols)  # +1: the heatmap panel
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 3.5 * n_rows))
    axes = np.atleast_1d(axes).ravel()

    for ax, field in zip(axes, fields):
        values = pd.to_numeric(df[field], errors="coerce").dropna()
        ylabel = "Count"
        xlabel = field
        if field in LOG_FIELDS:
            # Log X axis: log-spaced bins over the absolute values. Zeros
            # cannot be shown on a log axis and are dropped.
            if (values < 0).any():
                xlabel = f"|{field}|"
            values = values.abs()
            n_zeros = int((values == 0).sum())
            values = values[values > 0]
            if n_zeros > 0:
                xlabel += f" ({n_zeros} zeros not shown)"
            if len(values) > 0:
                bins = np.logspace(np.log10(values.min()),
                                   np.log10(values.max()), n_bins + 1)
                ax.hist(values, bins=bins, color="steelblue",
                        edgecolor="black", linewidth=0.3)
            ax.set_xscale("log")
        elif len(values) > 0:
            ax.hist(values, bins=n_bins, color="steelblue", edgecolor="black",
                    linewidth=0.3)
        ax.set_title(f"{field} (n={len(values)})")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    # 2D heatmap (linear scales, signed values) in the next free panel.
    n_used = len(fields)
    if n_used < len(axes):
        ax = axes[n_used]
        n_used += 1
        pair = df[["StartDistance", "PctSeqSimilarity"]].apply(
            pd.to_numeric, errors="coerce").dropna()
        # Log X axis: absolute StartDistance; zeros cannot be shown and are
        # dropped.
        pair["StartDistance"] = pair["StartDistance"].abs()
        n_zeros = int((pair["StartDistance"] == 0).sum())
        pair = pair[pair["StartDistance"] > 0]
        if len(pair) > 0:
            x = pair["StartDistance"]
            xbins = np.logspace(np.log10(x.min()), np.log10(x.max()),
                                n_bins + 1)
            h = ax.hist2d(x, pair["PctSeqSimilarity"],
                          bins=[xbins, n_bins], cmap="viridis",
                          norm=matplotlib.colors.LogNorm(), cmin=1)
            fig.colorbar(h[3], ax=ax, label="Count (log scale)")
        ax.set_xscale("log")
        ax.set_box_aspect(1)  # square axes
        ax.set_title(f"StartDistance vs PctSeqSimilarity (n={len(pair)})")
        ax.set_xlabel("|StartDistance|" +
                      (f" ({n_zeros} zeros not shown)" if n_zeros else ""))
        ax.set_ylabel("PctSeqSimilarity")
    for ax in axes[n_used:]:
        ax.axis("off")

    fig.suptitle(f"{in_tsv} ({len(df)} records)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
