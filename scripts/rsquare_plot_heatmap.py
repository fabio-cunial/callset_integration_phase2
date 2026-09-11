#!/usr/bin/env python3
"""Render a 2D density heatmap of R^2 (y) versus allele frequency (x) from the
BED emitted by `Rsquare.java`:

    CHROM, P, P+1, ID, type, length, AC, N, R^2, ACs, Ns, IDs

where AC/N describe the long variant and ACs/Ns/IDs describe the best-matching
short. This plot uses the AC and R^2 fields; type and length are read only to
back the --type / --min-length / --max-length filters.

The x axis is AF = AC / (2T), T being the total number of samples in the cohort
(the same definition `Rsquare.java` uses for its --min-af cutoff, where
minAC = ceil(2*totalNSamples*MIN_AF)). AF is extremely right-skewed, so the x
bins are log-spaced; R^2 is bounded in [0,1], so the y bins are equally spaced.

Binning happens in AC space, not AF space: AC is an integer, so the log-spaced
edges are rounded to integers and de-duplicated, and only then rescaled by 2T
for drawing. Every bin therefore covers a whole number of AC values, which
removes the empty/doubled stripes you otherwise get at low AF where a log bin
would be narrower than one allele.

Records with R^2 = -1 (the `Rsquare.java` sentinel for "no short variant to
compare against") are dropped and reported, never plotted as 0.

Usage:
    python plot_rsquare_heatmap.py rsquare.bed 6300
    python plot_rsquare_heatmap.py rsquare.bed 6300 --type DEL,INS --out di.png
    python plot_rsquare_heatmap.py rsquare.bed 6300 --x-bins 40 --dark
    python plot_rsquare_heatmap.py rsquare.bed 6300 --dump-counts  # table view

Colors: matplotlib's default sequential colormap, viridis (dark = near zero,
light = dense). Empty bins are drawn in the surface color, not in the colormap's
low step, so "no variants" stays distinct from "few variants".
"""

import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless: write files, never open a window
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import FuncFormatter, LogLocator

# Column indexes in the `Rsquare.java` output BED. CHROM (0), ID (3) and IDs
# (11) are strings and are never parsed.
BED_TYPE, BED_LENGTH, BED_AC, BED_N, BED_R2 = 4, 5, 6, 7, 8
BED_N_COLUMNS = 12
USECOLS = (BED_TYPE, BED_LENGTH, BED_AC, BED_N, BED_R2)

# Column indexes in the array load() returns (usecols compacts them).
COL_TYPE, COL_LENGTH, COL_AC, COL_N, COL_R2 = range(5)

# `Rsquare.java`: TYPE_SNP=0, TYPE_DEL=1, TYPE_INS=2, TYPE_SUB=3.
TYPE_NAMES = {0: "SNP", 1: "DEL", 2: "INS", 3: "SUB"}
TYPE_CODES = {v: k for k, v in TYPE_NAMES.items()}

# Axes-coordinate y of the x axis label. The legend is pinned to the same line.
XLABEL_Y = -0.12

THEMES = {
    "light": {"surface": "#fcfcfb", "primary": "#0b0b0b", "secondary": "#52514e"},
    "dark": {"surface": "#1a1a19", "primary": "#ffffff", "secondary": "#c3c2b7"},
}


def data_lines(path):
    """Yield the BED's data lines, skipping `#` comments, blank lines and an
    optional header row (a line whose R^2 field does not parse as a float).
    """
    checked = False
    with open(path) as f:
        for n, line in enumerate(f, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < BED_N_COLUMNS:
                sys.exit(f"error: {path}:{n} has {len(fields)} columns, expected "
                         f"{BED_N_COLUMNS} (CHROM,P,P+1,ID,type,length,AC,N,R^2,"
                         f"ACs,Ns,IDs)")
            if not checked:
                checked = True
                try:
                    float(fields[BED_R2])
                except ValueError:
                    continue  # header row
            yield line


def load(path):
    """Read the numeric columns of the `Rsquare.java` BED into an (n,5) array."""
    data = np.loadtxt(data_lines(path), delimiter="\t", usecols=USECOLS, ndmin=2)
    if data.size == 0:
        sys.exit(f"error: no data rows in {path}")
    return data


def log_integer_edges(lo, hi, n_bins):
    """Log-spaced bin edges over the integers [lo,hi], snapped to integers and
    de-duplicated so no bin is narrower than one AC value.

    Returns half-open edges e with e[0]=lo and e[-1]=hi+1, i.e. bin k holds the
    ACs in [e[k], e[k+1]-1].
    """
    raw = np.geomspace(lo, hi + 1, n_bins + 1)
    edges = np.unique(np.round(raw).astype(np.int64))
    edges[0], edges[-1] = lo, hi + 1
    return np.unique(edges)


def af_tick(v, _):
    """AF tick label: plain decimal where that is readable, 10^-k below it."""
    if v <= 0:
        return ""
    if v >= 0.01:
        return f"{v:g}"
    return f"$10^{{{int(round(np.log10(v)))}}}$"


def parse_types(spec):
    """`--type DEL,INS` or `--type 1,2` -> set of integer type codes."""
    out = set()
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok.upper() in TYPE_CODES:
            out.add(TYPE_CODES[tok.upper()])
        else:
            try:
                out.add(int(tok))
            except ValueError:
                sys.exit(f"error: unknown --type value {tok!r}")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bed", help="output BED of Rsquare.java: "
                                "CHROM,P,P+1,ID,type,length,AC,N,R^2,ACs,Ns,IDs")
    ap.add_argument("total_samples", type=int,
                    help="T, total number of samples in the cohort; the x axis "
                         "is AF = AC/(2T)")
    ap.add_argument("--out", help="output PNG (default: <bed>.af_r2_heatmap.png)")
    ap.add_argument("--x-bins", type=int, default=32, help="target number of log-spaced AC bins (default 32)")
    ap.add_argument("--y-bins", type=int, default=20, help="number of equal R^2 bins (default 20)")
    ap.add_argument("--max-ac", type=int, help="upper AC limit (default: max in the data)")
    ap.add_argument("--type", help="keep only these types, e.g. DEL,INS")
    ap.add_argument("--min-length", type=float, help="keep only variants of length >= this")
    ap.add_argument("--max-length", type=float, help="keep only variants of length <= this")
    ap.add_argument("--color-scale", choices=("auto", "log", "linear"), default="auto",
                    help="count scale for the color ramp (default auto: log when counts span >2 orders of magnitude)")
    ap.add_argument("--no-median", action="store_true", help="omit the per-AC-bin median R^2 trend line")
    ap.add_argument("--dark", action="store_true", help="render on the dark surface")
    ap.add_argument("--dump-counts", action="store_true",
                    help="also write the bin counts as a TSV table next to the PNG")
    ap.add_argument("--title", help="override the chart title")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()
    if args.total_samples < 1:
        sys.exit("error: total_samples must be >= 1")
    n_alleles = 2 * args.total_samples  # AF = AC / (2T)

    data = load(args.bed)
    n_total = len(data)

    # A wrong T (or a non-diploid callset) shows up here rather than as a
    # silently mislabelled axis.
    observed_max_ac = int(data[:, COL_AC].max())
    if observed_max_ac > n_alleles:
        sys.exit(f"error: {args.bed} has AC up to {observed_max_ac:,}, which exceeds "
                 f"2T = {n_alleles:,} for T = {args.total_samples:,} samples")

    # ---- filter ------------------------------------------------------------
    keep = np.ones(n_total, dtype=bool)
    if args.type:
        wanted = parse_types(args.type)
        keep &= np.isin(data[:, COL_TYPE].astype(np.int64), sorted(wanted))
    if args.min_length is not None:
        keep &= data[:, COL_LENGTH] >= args.min_length
    if args.max_length is not None:
        keep &= data[:, COL_LENGTH] <= args.max_length
    n_filtered_out = int((~keep).sum())

    uncomparable = keep & (data[:, COL_R2] < 0)  # R^2 == -1 sentinel
    n_uncomparable = int(uncomparable.sum())
    keep &= ~uncomparable
    keep &= data[:, COL_AC] >= 1

    ac = data[keep, COL_AC].astype(np.int64)
    r2 = np.clip(data[keep, COL_R2], 0.0, 1.0)
    if ac.size == 0:
        sys.exit("error: no records left to plot after filtering")

    max_ac = args.max_ac if args.max_ac else int(ac.max())
    in_range = ac <= max_ac
    n_over_max = int((~in_range).sum())
    ac, r2 = ac[in_range], r2[in_range]

    # ---- bin ---------------------------------------------------------------
    # Bin from the smallest AC actually present, not from 1: Rsquare.java's
    # --min-af cutoff leaves nothing below ceil(2T*min_af), and starting at 1
    # would spend most of a log axis on empty bins.
    min_ac = max(1, int(ac.min()))
    x_edges = log_integer_edges(min_ac, max_ac, args.x_bins)
    y_edges = np.linspace(0.0, 1.0, args.y_bins + 1)
    counts, _, _ = np.histogram2d(ac, r2, bins=[x_edges, y_edges])

    theme = THEMES["dark" if args.dark else "light"]
    # matplotlib's default sequential colormap. .copy() so set_bad() does not
    # mutate the globally registered instance.
    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad(theme["surface"])  # empty cells take the surface color
    grid = np.ma.masked_where(counts.T == 0, counts.T)

    vmax = counts.max()
    use_log = args.color_scale == "log" or (args.color_scale == "auto" and vmax >= 100)
    norm = LogNorm(vmin=1, vmax=max(vmax, 2)) if use_log else Normalize(vmin=0, vmax=vmax)

    # ---- draw --------------------------------------------------------------
    plt.rcParams.update({
        "font.size": 10,
        "text.color": theme["primary"],
        "axes.labelcolor": theme["secondary"],
        "xtick.color": theme["secondary"],
        "ytick.color": theme["secondary"],
    })
    fig, ax = plt.subplots(figsize=(9.0, 5.5), facecolor=theme["surface"])
    ax.set_facecolor(theme["surface"])

    # Bins are built on integer AC (above); rescale to AF only for drawing, so
    # integer AC k occupies [(k-0.5)/2T, (k+0.5)/2T).
    x_draw = (x_edges - 0.5) / n_alleles
    # A 1px surface gap between cells only when the grid is coarse enough for it
    # to read as separation rather than as noise.
    fine = (len(x_edges) - 1) * args.y_bins > 400
    mesh = ax.pcolormesh(x_draw, y_edges, grid, cmap=cmap, norm=norm,
                         edgecolors="none" if fine else theme["surface"],
                         linewidth=0 if fine else 0.5)

    if not args.no_median:
        centers, medians = [], []
        idx = np.digitize(ac, x_edges) - 1
        for k in range(len(x_edges) - 1):
            sel = idx == k
            if sel.sum() >= 10:  # don't draw a "median" of a handful of points
                centers.append(np.sqrt(max(x_edges[k], 0.5) * (x_edges[k + 1] - 0.5))
                               / n_alleles)
                medians.append(np.median(r2[sel]))
        if centers:
            # Thin dashed red: red is absent from viridis, so the line stays
            # readable over any cell without needing a halo to lift it off.
            ax.plot(centers, medians, color="red", linewidth=1.0, linestyle="--",
                    zorder=4, label="median")
            # Outside the axes, on the same line as the x label and flush with
            # the right spine, so it never covers a data cell.
            leg = ax.legend(loc="upper right", bbox_to_anchor=(1.0, XLABEL_Y),
                            bbox_transform=ax.transAxes, frameon=False,
                            handlelength=1.8, borderaxespad=0.0, borderpad=0.0)
            for text in leg.get_texts():
                text.set_color(theme["secondary"])

    ax.set_xscale("log")
    ax.set_xlim(x_draw[0], x_draw[-1])
    ax.set_ylim(0.0, 1.0)
    # A min-AF cutoff can leave the axis spanning barely a decade, where powers
    # of ten alone give one or two labels; add the 2/5 steps in that case.
    span_decades = np.log10(x_draw[-1] / x_draw[0])
    ax.xaxis.set_major_locator(
        LogLocator(base=10.0, subs=(1.0,) if span_decades >= 3 else (1.0, 2.0, 5.0)))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(np.arange(2, 10) * 0.1)))
    ax.xaxis.set_major_formatter(FuncFormatter(af_tick))
    ax.xaxis.set_minor_formatter(lambda v, _: "")
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_xlabel(f"AF ({args.total_samples:,} "
                  f"samples)")
    ax.xaxis.set_label_coords(0.5, XLABEL_Y)  # pins the label the legend aligns to
    ax.set_ylabel("max $R^2$ with a short call")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("bottom", "left"):
        ax.spines[side].set_color(theme["secondary"])
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(length=3, width=0.8)

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label("long calls per bin" + (" (log scale)" if use_log else ""),
                   color=theme["secondary"])
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(color=theme["secondary"], labelcolor=theme["secondary"], length=3, width=0.8)
    if use_log:
        cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))

    title = args.title or f"{os.path.basename(args.bed)}"
    ax.set_title(title, color=theme["primary"], loc="left", fontsize=12, pad=30)
    subtitle = f"{len(ac):,} total long calls"
    if n_filtered_out:
        subtitle += f" · {n_filtered_out:,} excluded by filters"
    if n_uncomparable:
        subtitle += f" · {n_uncomparable:,} with $R^2$=-1 (no comparable short) dropped"
    if n_over_max:
        subtitle += (f" · {n_over_max:,} above AF {max_ac / n_alleles:.3g} "
                     f"(AC {max_ac:,}) dropped")
    ax.text(0.0, 1.012, subtitle, transform=ax.transAxes, color=theme["secondary"],
            fontsize=9, va="bottom")

    out = args.out or f"{os.path.splitext(args.bed)[0]}.af_r2_heatmap.png"
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", facecolor=theme["surface"])
    plt.close(fig)
    print(f"wrote {out}  ({len(ac):,} variants, {len(x_edges)-1} AF bins x {args.y_bins} R^2 bins)")
    # Printing the range makes a T that does not match the BED obvious: with
    # Rsquare.java's --min-af cutoff, min AF should land on it.
    print(f"  AC {min_ac:,}..{max_ac:,}  ->  AF {min_ac/n_alleles:.5g}..{max_ac/n_alleles:.5g}"
          f"  (T = {args.total_samples:,}, 2T = {n_alleles:,})")

    if args.dump_counts:
        # Table view of exactly what the heatmap encodes. Bins are half-open in
        # AC, so both the integer AC range and the AF range are given.
        table = f"{os.path.splitext(out)[0]}.counts.tsv"
        with open(table, "w") as f:
            f.write("ac_lo\tac_hi\taf_lo\taf_hi\t" + "\t".join(
                f"r2_{y_edges[j]:.3f}_{y_edges[j+1]:.3f}" for j in range(args.y_bins)) + "\n")
            for k in range(len(x_edges) - 1):
                f.write(f"{x_edges[k]}\t{x_edges[k+1]-1}\t"
                        f"{x_edges[k]/n_alleles:.6g}\t{(x_edges[k+1]-1)/n_alleles:.6g}\t" +
                        "\t".join(str(int(c)) for c in counts[k]) + "\n")
        print(f"wrote {table}")


if __name__ == "__main__":
    main()
