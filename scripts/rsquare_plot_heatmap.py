#!/usr/bin/env python3
"""Render a figure from the BED emitted by `Rsquare.java`: a 3x3 grid of panels,
one row per stratum of the callset (all types, then DEL only, then INS only),
each row holding the same three panels. On the left of a row, a 2D density
heatmap of R^2 (y) versus allele frequency (x); in the middle, the distribution
of all the R^2 values plotted in that row; on the right, the same kind of
heatmap with call length (x) in place of allele frequency. The BED columns are:

    CHROM, P, P+1, ID, type, length, AC, N, R^2, ACs, Ns, IDs

where AC/N describe the long variant and ACs/Ns/IDs describe the best-matching
short. This plot uses the type, length, AC and R^2 fields.

Within a row, all three panels describe exactly the same set of calls and share
the y axis and the R^2 bins, so the middle panel is the row marginal of either
heatmap (the same calls, summed over AF or over length); it answers "how are
R^2 values distributed overall", which a heatmap hides by spreading each row
over its x axis and by compressing counts through a color ramp.

Across rows, everything is shared as well: the same bin edges, the same axis
limits, one color scale and one colorbar. The DEL and INS rows are subsets of
the first row, so a cell, a bar or a color means the same number of calls
wherever it appears, and the three rows can be read as a decomposition of the
first rather than as three unrelated figures. Empty bins stay visible as
surface, which is what makes a stratum's missing regions legible.

The x axis of the left heatmaps is AF = AC / (2T), T the total number of samples
(the same definition `Rsquare.java` uses for its --min-af cutoff, where
minAC = ceil(2*totalNSamples*MIN_AF)); the x axis of the right heatmaps is the
length in bp that `Rsquare.java` computed for the long call. Both are extremely
right-skewed, so both are log-spaced; R^2 is bounded in [0,1], so the y bins are
equally spaced.

Binning happens in integer space, not in AF space: AC and length are integers,
so the log-spaced edges are rounded to integers and de-duplicated, and only the
AC edges are then rescaled by 2T for drawing. Every bin therefore covers a whole
number of AC (resp. length) values, which removes the empty/doubled stripes you
otherwise get at low AF where a log bin would be narrower than one allele.

Records with R^2 = -1 (the `Rsquare.java` sentinel for "no short variant to
compare against") are dropped and reported, never plotted as 0.

Usage:
    python plot_rsquare_heatmap.py rsquare.bed 6300
    python plot_rsquare_heatmap.py rsquare.bed 6300 --type DEL,INS --out di.png
    python plot_rsquare_heatmap.py rsquare.bed 6300 --x-bins 40 --dark
    python plot_rsquare_heatmap.py rsquare.bed 6300 --dump-counts  # table view

Colors: the heatmaps use matplotlib's default sequential colormap, viridis (dark
= near zero, light = dense). Empty bins are drawn in the surface color, not in
the colormap's low step, so "no variants" stays distinct from "few variants".
The middle panels are a single series, so they take one flat blue.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless: write files, never open a window
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator

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

# One row of panels per entry, top to bottom: (label, type code to keep). A
# `None` code keeps whatever survived the command-line filters, so the first row
# is the whole plotted callset and the others decompose it.
STRATA = (("All types", None), ("DEL", TYPE_CODES["DEL"]), ("INS", TYPE_CODES["INS"]))

# Axes-coordinate y of the x axis label. The legends are pinned to the same line.
XLABEL_Y = -0.16

# Axes-coordinate y of the per-panel header line, and of the figure-wide
# subtitle that sits one line above it on the top-left panel.
HEADER_Y = 1.02
SUBTITLE_Y = 1.10

Y_LABEL = "max $R^2$ with a short call"

THEMES = {
    # "bar" is the single-series color of the middle panels, stepped per surface.
    "light": {"surface": "#fcfcfb", "primary": "#0b0b0b", "secondary": "#52514e",
              "bar": "#2a78d6"},
    "dark": {"surface": "#1a1a19", "primary": "#ffffff", "secondary": "#c3c2b7",
             "bar": "#3987e5"},
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
    de-duplicated so no bin is narrower than one unit.

    Returns half-open edges e with e[0]=lo and e[-1]=hi+1, i.e. bin k holds the
    values in [e[k], e[k+1]-1].
    """
    raw = np.geomspace(lo, hi + 1, n_bins + 1)
    edges = np.unique(np.round(raw).astype(np.int64))
    edges[0], edges[-1] = lo, hi + 1
    return np.unique(edges)


def af_tick(v, _):
    """AF tick label: plain decimal where that is readable, m x 10^-k below it.

    The locator also places ticks at the 2 and 5 steps when the axis spans less
    than three decades, so the mantissa has to survive into the label: rounding
    log10 would print both 0.002 and 0.005 as a power of ten, and at the wrong
    decade.
    """
    if v <= 0:
        return ""
    if v >= 0.01:
        return f"{v:g}"
    exponent = int(np.floor(np.log10(v) + 1e-9))
    mantissa = v / 10.0 ** exponent
    if abs(mantissa - 1.0) < 1e-6:
        return f"$10^{{{exponent}}}$"
    return f"${mantissa:.3g}\\times10^{{{exponent}}}$"


def count_tick(v, _):
    """Count tick label, abbreviated: 1200 -> 1.2k, 3400000 -> 3.4M. Also used
    for the bp ticks of the length axis, which needs the same abbreviation.
    """
    if v < 0:
        return ""
    for scale, suffix in ((1e6, "M"), (1e3, "k")):
        if v >= scale:
            return f"{v / scale:g}{suffix}"
    return f"{v:g}"


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


def median_trend(values, r2, edges, scale=1.0, min_count=10):
    """Per-x-bin median R^2, as (centers, medians) ready to plot.

    `values` are the integers the bins were built on and `edges` their half-open
    integer edges; `scale` divides the centers into the axis' draw coordinates
    (2T for the AF axis, 1 for the length axis). A center is the geometric mid
    of the bin, to sit where the eye puts it on a log axis. Bins holding fewer
    than `min_count` calls are skipped: a "median" of a handful of points is
    noise drawn at the same weight as the rest of the line.
    """
    centers, medians = [], []
    idx = np.digitize(values, edges) - 1
    for k in range(len(edges) - 1):
        sel = idx == k
        if sel.sum() >= min_count:
            centers.append(np.sqrt(max(edges[k], 0.5) * (edges[k + 1] - 0.5)) / scale)
            medians.append(np.median(r2[sel]))
    return centers, medians


def style_spines(ax, theme, visible=("bottom", "left")):
    """Drop the top/right spines and tone the remaining ones down to the
    surface's secondary color."""
    for side in ("top", "right", "bottom", "left"):
        if side in visible:
            ax.spines[side].set_color(theme["secondary"])
            ax.spines[side].set_linewidth(0.8)
        else:
            ax.spines[side].set_visible(False)


def draw_heatmap(ax, x_draw, y_edges, counts, cmap, norm, theme):
    """Draw one R^2-versus-x density heatmap and return its QuadMesh.

    `counts` is the (x,y) histogram; it is transposed here because pcolormesh
    wants (row=y, col=x). A 1px surface gap between cells is drawn only when the
    grid is coarse enough for it to read as separation rather than as noise.
    """
    grid = np.ma.masked_where(counts.T == 0, counts.T)
    fine = counts.size > 400
    return ax.pcolormesh(x_draw, y_edges, grid, cmap=cmap, norm=norm,
                         edgecolors="none" if fine else theme["surface"],
                         linewidth=0 if fine else 0.5)


def style_log_x(ax, x_draw, formatter):
    """Log x axis clamped to the drawn bins. A min-AF cutoff (or a narrow length
    range) can leave the axis spanning barely a decade, where powers of ten
    alone give one or two labels; add the 2/5 steps in that case."""
    ax.set_xscale("log")
    ax.set_xlim(x_draw[0], x_draw[-1])
    span_decades = np.log10(x_draw[-1] / x_draw[0])
    ax.xaxis.set_major_locator(
        LogLocator(base=10.0, subs=(1.0,) if span_decades >= 3 else (1.0, 2.0, 5.0)))
    ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=tuple(np.arange(2, 10) * 0.1)))
    ax.xaxis.set_major_formatter(FuncFormatter(formatter))
    ax.xaxis.set_minor_formatter(lambda v, _: "")


def row_header(ax, theme, name, detail):
    """Header line of a row of panels: the stratum's name in bold, followed by
    `detail` in the regular weight, so the name alone carries the emphasis.

    The two runs are packed rather than placed at fixed x coordinates: the name
    is one to four words long depending on --type, and only the packer knows how
    wide it renders.
    """
    props = dict(color=theme["primary"], fontsize=10)
    box = HPacker(pad=0, sep=4, align="baseline", children=[
        TextArea(name, textprops=dict(props, fontweight="bold")),
        TextArea(detail, textprops=props)])
    ax.add_artist(AnchoredOffsetbox(loc="lower left", child=box, pad=0.0,
                                    borderpad=0.0, frameon=False,
                                    bbox_to_anchor=(0.0, HEADER_Y),
                                    bbox_transform=ax.transAxes))


def median_legend(ax, theme):
    """Legend for the median line, placed outside the axes on the same line as
    the x label and flush with the right spine, so it never covers a data cell.
    """
    leg = ax.legend(loc="upper right", bbox_to_anchor=(1.0, XLABEL_Y),
                    bbox_transform=ax.transAxes, frameon=False,
                    handlelength=1.8, borderaxespad=0.0, borderpad=0.0)
    for text in leg.get_texts():
        text.set_color(theme["secondary"])


def dump_counts(path, counts, x_edges, y_edges, x_lo_name, x_hi_name, extra=None):
    """Table view of exactly what one heatmap encodes: one row per x bin, then a
    final `all` row holding the column totals (the middle panel of that stratum).
    Bins are half-open, so the printed high edge is the last value the bin
    contains.

    `extra` maps an extra column name to a callable taking the bin index (or
    `None` for the totals row), for the AF columns the AC table also carries.
    """
    extra = extra or {}
    with open(path, "w") as f:
        f.write(f"{x_lo_name}\t{x_hi_name}\t" + "".join(f"{k}\t" for k in extra) +
                "\t".join(f"r2_{y_edges[j]:.3f}_{y_edges[j+1]:.3f}"
                          for j in range(len(y_edges) - 1)) + "\n")
        for k in range(len(x_edges) - 1):
            f.write(f"{x_edges[k]}\t{x_edges[k+1]-1}\t" +
                    "".join(f"{fn(k)}\t" for fn in extra.values()) +
                    "\t".join(str(int(c)) for c in counts[k]) + "\n")
        f.write("all\tall\t" + "".join(f"{fn(None)}\t" for fn in extra.values()) +
                "\t".join(str(int(c)) for c in counts.sum(axis=0)) + "\n")
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bed", help="output BED of Rsquare.java: "
                                "CHROM,P,P+1,ID,type,length,AC,N,R^2,ACs,Ns,IDs")
    ap.add_argument("total_samples", type=int,
                    help="T, total number of samples in the cohort; the x axis "
                         "of the left heatmaps is AF = AC/(2T)")
    ap.add_argument("--out", help="output PNG (default: <bed>.af_r2_heatmap.png)")
    ap.add_argument("--x-bins", type=int, default=32,
                    help="target number of log-spaced bins on the x axis of "
                         "each heatmap, i.e. AC bins on the left and length "
                         "bins on the right (default 32)")
    ap.add_argument("--y-bins", type=int, default=20,
                    help="number of equal R^2 bins, in all panels (default 20)")
    ap.add_argument("--max-ac", type=int, help="upper AC limit (default: max in the data)")
    ap.add_argument("--type", help="keep only these types, e.g. DEL,INS; this "
                                   "restricts the first row of panels, the DEL "
                                   "and INS rows being what they say regardless")
    ap.add_argument("--min-length", type=float, help="keep only variants of length >= this")
    ap.add_argument("--max-length", type=float, help="keep only variants of length <= this")
    ap.add_argument("--color-scale", choices=("auto", "log", "linear"), default="auto",
                    help="count scale for the color ramp shared by every "
                         "heatmap (default auto: log when counts span >2 "
                         "orders of magnitude)")
    ap.add_argument("--no-median", action="store_true",
                    help="omit the per-bin median R^2 trend line of the "
                         "heatmaps (and the overall median in the middle panels)")
    ap.add_argument("--dark", action="store_true", help="render on the dark surface")
    ap.add_argument("--dump-counts", action="store_true",
                    help="also write the bin counts of every heatmap as TSV "
                         "tables next to the PNG")
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
    # A length of 0 cannot be placed on a log axis. `Rsquare.java` never emits
    # one (getLength is >=1 for every type it assigns), so this only guards a
    # hand-edited BED, but dropping the record from *every* panel is what keeps
    # the middle panels the exact marginals of their heatmaps.
    n_nonpositive_length = int((keep & (data[:, COL_LENGTH] < 1)).sum())
    keep &= data[:, COL_LENGTH] >= 1

    types = data[keep, COL_TYPE].astype(np.int64)
    ac = data[keep, COL_AC].astype(np.int64)
    length = data[keep, COL_LENGTH].astype(np.int64)
    r2 = np.clip(data[keep, COL_R2], 0.0, 1.0)
    if ac.size == 0:
        sys.exit("error: no records left to plot after filtering")

    max_ac = args.max_ac if args.max_ac else int(ac.max())
    in_range = ac <= max_ac
    n_over_max = int((~in_range).sum())
    types, ac, length, r2 = types[in_range], ac[in_range], length[in_range], r2[in_range]

    # ---- bin ---------------------------------------------------------------
    # Bin from the smallest value actually present, not from 1: Rsquare.java's
    # --min-af cutoff leaves nothing below ceil(2T*min_af), and its
    # --min-long-length cutoff leaves nothing below that length; starting at 1
    # would spend most of a log axis on empty bins.
    #
    # The edges come from the whole callset and are then reused for every
    # stratum, so the rows stack into one comparable grid: where a stratum has
    # no calls, its panel shows surface rather than a rebinned axis of its own.
    min_ac = max(1, int(ac.min()))
    x_edges = log_integer_edges(min_ac, max_ac, args.x_bins)
    y_edges = np.linspace(0.0, 1.0, args.y_bins + 1)
    min_length, max_length = max(1, int(length.min())), int(length.max())
    len_edges = log_integer_edges(min_length, max_length, args.x_bins)

    strata = []
    for label, code in STRATA:
        sel = np.ones(len(ac), dtype=bool) if code is None else types == code
        strata.append({
            "label": label if code is not None or not args.type
                     else f"all selected types ({args.type.upper()})",
            "ac": ac[sel], "length": length[sel], "r2": r2[sel],
            "counts": np.histogram2d(ac[sel], r2[sel], bins=[x_edges, y_edges])[0],
            "len_counts": np.histogram2d(length[sel], r2[sel],
                                         bins=[len_edges, y_edges])[0],
        })

    theme = THEMES["dark" if args.dark else "light"]
    # matplotlib's default sequential colormap. .copy() so set_bad() does not
    # mutate the globally registered instance.
    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad(theme["surface"])  # empty cells take the surface color

    # One norm for every heatmap, so that equally dark cells mean equally many
    # calls anywhere in the figure and a single colorbar can describe them all.
    vmax = max(max(s["counts"].max(), s["len_counts"].max()) for s in strata)
    use_log = args.color_scale == "log" or (args.color_scale == "auto" and vmax >= 100)
    norm = LogNorm(vmin=1, vmax=max(vmax, 2)) if use_log else Normalize(vmin=0, vmax=vmax)
    # Likewise for the middle column: the strata are subsets of the first row,
    # so a shared count axis makes the bars show composition and not just shape.
    max_bar = max(s["counts"].sum(axis=0).max() for s in strata)

    # ---- draw --------------------------------------------------------------
    plt.rcParams.update({
        "font.size": 10,
        "text.color": theme["primary"],
        "axes.labelcolor": theme["secondary"],
        "xtick.color": theme["secondary"],
        "ytick.color": theme["secondary"],
    })
    # One row of panels per stratum. Within a row: the AF heatmap and, glued to
    # it, its marginal; then the same heatmap over call length. The empty third
    # gridspec column is the gutter that separates the two halves; without it
    # the marginal would read as belonging to whichever heatmap it sits closer
    # to. Every panel shares the R^2 axis and every column shares its x axis, so
    # a band, a column or a cell lines up everywhere it appears.
    fig = plt.figure(figsize=(17.0, 4.6 * len(strata)), facecolor=theme["surface"])
    gs = fig.add_gridspec(len(strata), 4, width_ratios=[3.2, 1.0, 0.55, 3.2],
                          wspace=0.05, hspace=0.3)
    anchors = {}  # column -> the axes the rest of that column shares x with
    for s, stratum in enumerate(strata):
        panels = []
        for col in (0, 1, 3):
            ax = fig.add_subplot(gs[s, col], sharex=anchors.get(col),
                                 sharey=anchors.get(0))
            anchors.setdefault(col, ax)
            ax.set_facecolor(theme["surface"])
            panels.append(ax)
        stratum["axes"] = panels

    # AC bins are built on integers (above); rescale to AF only for drawing, so
    # integer AC k occupies [(k-0.5)/2T, (k+0.5)/2T). Length is already in the
    # units of its axis and only needs the same half-unit centering.
    x_draw = (x_edges - 0.5) / n_alleles
    len_draw = len_edges - 0.5

    # The figure's one median legend goes on the lowest row that has a line to
    # explain, which is the bottom row unless a stratum came out empty.
    legend_row = max((i for i, st in enumerate(strata) if len(st["ac"])), default=-1)

    mesh = None
    for s, stratum in enumerate(strata):
        ax, ax_dist, ax_len = stratum["axes"]
        bottom = s == len(strata) - 1
        n_calls = len(stratum["ac"])

        heat = draw_heatmap(ax, x_draw, y_edges, stratum["counts"], cmap, norm, theme)
        mesh = mesh or heat  # the norm is shared, so any mesh describes the bar
        draw_heatmap(ax_len, len_draw, y_edges, stratum["len_counts"], cmap, norm, theme)
        if n_calls == 0:
            # An empty stratum is a legitimate result (a BED with no INS, or a
            # --type that excludes one), and an unexplained blank panel is not.
            for panel in (ax, ax_len):
                panel.text(0.5, 0.5, "no calls", transform=panel.transAxes,
                           ha="center", va="center", color=theme["secondary"],
                           fontsize=10)

        if not args.no_median:
            for panel, values, edges, scale in (
                    (ax, stratum["ac"], x_edges, n_alleles),
                    (ax_len, stratum["length"], len_edges, 1.0)):
                centers, medians = median_trend(values, stratum["r2"], edges, scale)
                if centers:
                    # Thin dashed red: red is absent from viridis, so the line
                    # stays readable over any cell without needing a halo to
                    # lift it off.
                    panel.plot(centers, medians, color="red", linewidth=1.0,
                               linestyle="--", zorder=4, label="median per bin")
                    # One legend for the whole figure, on the line the bottom
                    # row's x labels sit on; the encoding is the same in every
                    # row. "per bin" needs no further qualification there, and
                    # stays short enough to clear the x label it shares a line
                    # with.
                    if s == legend_row:
                        median_legend(panel, theme)

        # ---- left panel: R^2 versus AF -------------------------------------
        style_log_x(ax, x_draw, af_tick)
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.set_ylabel(Y_LABEL)
        style_spines(ax, theme)
        ax.tick_params(length=3, width=0.8, labelbottom=bottom)
        if bottom:
            ax.set_xlabel(f"AF ({args.total_samples:,} samples)")
            ax.xaxis.set_label_coords(0.5, XLABEL_Y)  # the line the legend aligns to

        # ---- middle panel: the distribution of this row's R^2 values --------
        # Summing the heatmap's columns rather than re-histogramming r2 keeps
        # the panels of a row describing the exact same calls, by construction.
        r2_totals = stratum["counts"].sum(axis=0)
        ax_dist.barh(y_edges[:-1], r2_totals, height=np.diff(y_edges), align="edge",
                     color=theme["bar"], edgecolor=theme["surface"], linewidth=0.8,
                     zorder=2)
        if not args.no_median and n_calls:
            # The overall median, i.e. where this distribution's mass splits;
            # the heatmaps' lines are the same statistic taken per AF / per
            # length bin.
            median_all = float(np.median(stratum["r2"]))
            ax_dist.axhline(median_all, color="red", linewidth=1.0, linestyle="--",
                            zorder=4)
            # Flush right, just above the line; on a surface patch because the
            # bar it sits over is long whenever the median lands near the mode.
            ax_dist.text(0.97, median_all + 0.012, f"median {median_all:.2f}",
                         transform=ax_dist.get_yaxis_transform(), ha="right",
                         va="bottom", fontsize=9, color=theme["secondary"], zorder=5,
                         bbox=dict(facecolor=theme["surface"], edgecolor="none", pad=1.5))
        ax_dist.set_xlim(0, max(max_bar, 1) * 1.02)
        ax_dist.xaxis.set_major_locator(MaxNLocator(3, integer=True))
        ax_dist.xaxis.set_major_formatter(FuncFormatter(count_tick))
        style_spines(ax_dist, theme)  # the left spine is the bars' baseline
        ax_dist.tick_params(length=3, width=0.8, left=False, labelleft=False,
                            labelbottom=bottom)
        if bottom:
            ax_dist.set_xlabel("long calls")
            ax_dist.xaxis.set_label_coords(0.5, XLABEL_Y)  # same line as the AF label

        # ---- right panel: the same calls over call length -------------------
        style_log_x(ax_len, len_draw, count_tick)
        # The gutter puts this panel far enough from the left one that reading
        # its rows off those y tick labels no longer works, so it repeats them.
        ax_len.set_yticks(np.linspace(0, 1, 6))
        ax_len.set_ylabel(Y_LABEL)
        style_spines(ax_len, theme)
        ax_len.tick_params(length=3, width=0.8, labelbottom=bottom)
        if bottom:
            ax_len.set_xlabel("call length (bp)")
            ax_len.xaxis.set_label_coords(0.5, XLABEL_Y)

        # ---- headers --------------------------------------------------------
        # The stratum name sits over its leftmost panel; the three panels of a
        # row are aligned under it, so one label names the whole row.
        row_header(ax, theme, stratum["label"], f"· {n_calls:,} long calls")
        if s == 0:
            # Naming the other two panels once, on the top row, is enough: the
            # rows below repeat the same layout.
            ax_dist.text(0.0, HEADER_Y, "all $R^2$ values",
                         transform=ax_dist.transAxes, color=theme["secondary"],
                         fontsize=9, va="bottom")
            ax_len.text(0.0, HEADER_Y, "The same calls, by length.",
                        transform=ax_len.transAxes, color=theme["secondary"],
                        fontsize=9, va="bottom")

    # The colorbar spans the full height of the length column, since one norm
    # covers every heatmap in the figure.
    cbar = fig.colorbar(mesh, ax=[s["axes"][2] for s in strata], pad=0.04)
    cbar.set_label("long calls per bin" + (" (log scale)" if use_log else ""),
                   color=theme["secondary"])
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(color=theme["secondary"], labelcolor=theme["secondary"], length=3, width=0.8)
    if use_log:
        cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))

    top_ax = strata[0]["axes"][0]
    title = args.title or f"{os.path.basename(args.bed)}"
    top_ax.set_title(title, color=theme["primary"], loc="left", fontsize=12, pad=46)
    subtitle = f"{len(ac):,} long calls plotted"
    if n_filtered_out:
        subtitle += f" · {n_filtered_out:,} excluded by filters"
    if n_uncomparable:
        subtitle += f" · {n_uncomparable:,} with $R^2$=-1 (no comparable short) dropped"
    if n_nonpositive_length:
        subtitle += f" · {n_nonpositive_length:,} with length <1 dropped"
    if n_over_max:
        subtitle += (f" · {n_over_max:,} above AF {max_ac / n_alleles:.3g} "
                     f"(AC {max_ac:,}) dropped")
    top_ax.text(0.0, SUBTITLE_Y, subtitle, transform=top_ax.transAxes,
                color=theme["secondary"], fontsize=9, va="bottom")

    out = args.out or f"{os.path.splitext(args.bed)[0]}.af_r2_heatmap.png"
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", facecolor=theme["surface"])
    plt.close(fig)
    print(f"wrote {out}  ({len(ac):,} variants, {len(x_edges)-1} AF bins and "
          f"{len(len_edges)-1} length bins x {args.y_bins} R^2 bins)")
    print("  " + " · ".join(f"{s['label']}: {len(s['ac']):,}" for s in strata))
    # Printing the ranges makes a T that does not match the BED obvious: with
    # Rsquare.java's --min-af cutoff, min AF should land on it.
    print(f"  AC {min_ac:,}..{max_ac:,}  ->  AF {min_ac/n_alleles:.5g}..{max_ac/n_alleles:.5g}"
          f"  (T = {args.total_samples:,}, 2T = {n_alleles:,})")
    print(f"  length {min_length:,}..{max_length:,} bp")

    if args.dump_counts:
        # Two tables per stratum, one per heatmap; each carries a final `all`
        # row, which is that stratum's middle panel.
        stem = os.path.splitext(out)[0]
        for (_, code), stratum in zip(STRATA, strata):
            tag = "all" if code is None else TYPE_NAMES[code]
            dump_counts(f"{stem}.{tag}.counts.tsv", stratum["counts"], x_edges,
                        y_edges, "ac_lo", "ac_hi",
                        extra={"af_lo": lambda k: f"{x_edges[0 if k is None else k]/n_alleles:.6g}",
                               "af_hi": lambda k: f"{(x_edges[-1 if k is None else k+1]-1)/n_alleles:.6g}"})
            dump_counts(f"{stem}.{tag}.length_counts.tsv", stratum["len_counts"],
                        len_edges, y_edges, "length_lo", "length_hi")


if __name__ == "__main__":
    main()
