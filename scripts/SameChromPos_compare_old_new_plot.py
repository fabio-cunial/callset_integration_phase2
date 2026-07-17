#!/usr/bin/env python3
"""
Plot a 3-panel figure comparing an OLD and a NEW SV callset.

Panel 1: one circle per row of the counts CSV, showing three per-row ratios
         (OLD_ONLY fraction, NEW_ONLY fraction, NEW_ONLY_IN_TRS fraction) at
         three distinct X coordinates, with horizontal jitter and an overlaid
         box-and-whiskers plot.
Panel 2: histogram of SVLEN for SVTYPE=INS records.
Panel 3: histogram of SVLEN for SVTYPE=DEL records.

Usage:
    SameChromPos_compare_old_new_plot.py counts.csv svtype_svlen.csv title

The counts CSV has columns OLD_ONLY,NEW_ONLY,NEW_ONLY_IN_TRS,COMMON.
The svtype/svlen CSV has columns SVTYPE,SVLEN.
Both files may or may not carry a header line.
"""

import sys
import csv

import numpy as np
import matplotlib.pyplot as plt


COUNTS_COLS = ["OLD_ONLY", "NEW_ONLY", "NEW_ONLY_IN_TRS", "COMMON"]


def _is_number(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


def load_counts(path):
    """Return a dict of column_name -> np.array of floats."""
    with open(path, newline="") as f:
        rows = [r for r in csv.reader(f) if r and any(c.strip() for c in r)]
    if not rows:
        raise SystemExit(f"No data in {path}")
    # Skip a header line if the first row is not fully numeric.
    if not all(_is_number(c) for c in rows[0][:4]):
        rows = rows[1:]
    data = np.array([[float(c) for c in r[:4]] for r in rows], dtype=float)
    return {name: data[:, i] for i, name in enumerate(COUNTS_COLS)}


def load_svtype_svlen(path):
    """Return dict svtype -> np.array of |SVLEN|."""
    result = {}
    with open(path, newline="") as f:
        for r in csv.reader(f):
            if not r or len(r) < 2:
                continue
            svtype, svlen = r[0].strip(), r[1].strip()
            if not _is_number(svlen):
                continue  # header or missing value
            result.setdefault(svtype, []).append(abs(float(svlen)))
    return {k: np.array(v, dtype=float) for k, v in result.items()}


def safe_ratio(num, den):
    """Elementwise num/den, returning NaN where den == 0."""
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full(num.shape, np.nan)
    mask = den != 0
    out[mask] = num[mask] / den[mask]
    return out


def plot_ratios(ax, counts):
    old_only = counts["OLD_ONLY"]
    new_only = counts["NEW_ONLY"]
    new_only_trs = counts["NEW_ONLY_IN_TRS"]
    common = counts["COMMON"]

    ratios = [
        safe_ratio(old_only, old_only + common),
        safe_ratio(new_only, new_only + common),
        safe_ratio(new_only_trs, new_only),
    ]
    labels = [
        "Old only / Old",
        "New only / New",
        "New only in TRs / New only",
    ]
    positions = [1, 2, 3]

    # Box-and-whiskers underneath the points (drop NaNs per group).
    box_data = [r[~np.isnan(r)] for r in ratios]
    ax.boxplot(
        box_data,
        positions=positions,
        widths=0.5,
        showfliers=False,
        patch_artist=True,
        boxprops=dict(facecolor="#dddddd", alpha=0.6),
        medianprops=dict(color="black"),
        zorder=1,
    )

    rng = np.random.default_rng(0)
    for pos, r in zip(positions, ratios):
        valid = r[~np.isnan(r)]
        jitter = rng.uniform(-0.12, 0.12, size=valid.shape)
        ax.scatter(
            np.full(valid.shape, pos) + jitter,
            valid,
            s=40,
            facecolors="none",
            edgecolors="tab:blue",
            linewidths=1.2,
            zorder=2,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Fraction")
    ax.set_title("New vs Old pipeline, 15x, 10 HPRC/HGSVC samples")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(axis="y", linestyle=":", alpha=0.5)


def plot_svlen_hist(ax, values, svtype, color):
    if values.size == 0:
        ax.text(0.5, 0.5, f"No {svtype} records", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title(f"SVLEN histogram ({svtype})")
        return
    weights = np.full(values.shape, 1.0 / values.size)
    ax.hist(values, bins=100, weights=weights, color=color,
            edgecolor="black", linewidth=0.3)
    ax.set_xlabel("SVLEN (bp)")
    ax.set_ylabel("Fraction")
    ax.set_title(f"New only, {svtype}, n={values.size}")
    ax.grid(axis="y", linestyle=":", alpha=0.5)


def main():
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    counts_csv = sys.argv[1]
    svlen_csv = sys.argv[2]
    title = sys.argv[3]

    counts = load_counts(counts_csv)
    svlen_by_type = load_svtype_svlen(svlen_csv)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    plot_ratios(axes[0], counts)
    plot_svlen_hist(axes[1], svlen_by_type.get("INS", np.array([])), "INS", "tab:green")
    plot_svlen_hist(axes[2], svlen_by_type.get("DEL", np.array([])), "DEL", "tab:red")

    # Prepend the given title to every panel's title.
    for ax in axes:
        ax.set_title(f"{title}, {ax.get_title()}")

    fig.tight_layout()
    plt.show()



if __name__ == "__main__":
    main()
