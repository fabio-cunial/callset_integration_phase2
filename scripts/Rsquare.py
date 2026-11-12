#!/usr/bin/env python3
"""Render a figure from two BEDs emitted by `Rsquare.java`, one holding the long
calls inside tandem repeats and one holding those outside (see `Rsquare.sh`),
restricted to one SV type (DEL or INS) or not (ALL). The figure is a 2x2 grid of square
panels: the left column describes the calls inside TRs, the right column those
outside TRs. The top row is a 2D density heatmap of R^2 (y) versus allele
frequency (x), with a narrow marginal histogram of all the R^2 values of the
panel glued to its right; the bottom row is the same kind of heatmap with call
length (x) in place of allele frequency. The BED columns are:

    CHROM, P, P+1, ID, type, length, AC, N, R^2, ACs, Ns, IDs

where AC/N describe the long variant and ACs/Ns/IDs describe the best-matching
short. This plot uses the type, length, AC and R^2 fields.

Within a column, both heatmaps and the marginal describe exactly the same set of
calls and share the y axis and the R^2 bins, so the marginal is the row marginal
of either heatmap (the same calls, summed over AF or over length); it answers
"how are R^2 values distributed overall", which a heatmap hides by spreading
each row over its x axis and by compressing counts through a color ramp. The
marginal shows the fraction of the column's calls in each R^2 bin, so it sums
to 1; it has no axis and is scaled to its own largest bar, since it is read
only for its shape.

Across columns, everything is shared as well: the same bin edges, the same axis
limits, one color scale. A cell or a color therefore means the same number of
calls in both columns, and the two can be compared directly. Empty bins stay visible as surface, which is what makes a stratum's
missing regions legible.

The x axis of the top heatmaps is AF = AC / (2T), T the total number of samples
(the same definition `Rsquare.java` uses for its --min-af cutoff, where
minAC = ceil(2*totalNSamples*MIN_AF)); the x axis of the bottom heatmaps is the
length in bp that `Rsquare.java` computed for the long call. Both are extremely
right-skewed, so both are log-spaced; R^2 is bounded in [0,1], so the y bins are
equally spaced.

Binning happens in integer space, not in AF space: AC and length are integers,
so the log-spaced edges are rounded to integers and de-duplicated, and only the
AC edges are then rescaled by 2T for drawing. Every bin therefore covers a whole
number of AC (resp. length) values, which removes the empty/doubled stripes you
otherwise get at low AF where a log bin would be narrower than one allele.

Records with R^2 = -1 (the `Rsquare.java` sentinel for "no short variant to
compare against") are dropped and reported on stdout, never plotted as 0.

Usage:
    python Rsquare.py inside_tr.bed outside_tr.bed 6300 DEL
    python Rsquare.py inside_tr.bed outside_tr.bed 6300 all
    python Rsquare.py inside_tr.bed outside_tr.bed 6300 INS --out ins.png
    python Rsquare.py inside_tr.bed outside_tr.bed 6300 DEL --title-prefix "DEL,"
    python Rsquare.py inside_tr.bed outside_tr.bed 6300 DEL --x-bins 40 --dark
    python Rsquare.py inside_tr.bed outside_tr.bed 6300 DEL --dump-counts

Colors: the heatmaps use matplotlib's default sequential colormap, viridis (dark
= near zero, light = dense). Empty bins are drawn in the surface color, not in
the colormap's low step, so "no variants" stays distinct from "few variants".
The marginals are a single series, so they take one flat blue. Medians are
drawn as a continuous white line.
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
SVTYPE_ALL = "ALL"  # svtype that keeps every record, whatever its type
SVTYPES = ("DEL", "INS", SVTYPE_ALL)

# Axes-coordinate y of the per-column header line.
HEADER_Y = 1.02

# Figure geometry, in inches. Laying the axes out by hand rather than through a
# gridspec is what keeps every heatmap exactly square while its marginal stays
# glued to it at a fixed width.
PANEL = 3.4          # side of a square heatmap
MARGINAL = 0.4       # width of a marginal histogram
MARGINAL_GAP = 0.04  # between a heatmap and its marginal
GUTTER = 0.75        # between the two columns; holds the right y tick labels
ROW_GAP = 0.85       # between the two rows; holds the top x tick labels/label
MARGIN_LEFT, MARGIN_RIGHT, MARGIN_BOTTOM, MARGIN_TOP = 0.8, 0.1, 0.7, 0.4

MEDIAN_COLOR = "white"

Y_LABEL = "$R^2$"

THEMES = {
    # "bar" is the single-series color of the marginals, stepped per surface.
    # "border" frames every panel; black would vanish on the dark surface.
    "light": {"surface": "#fcfcfb", "primary": "#0b0b0b", "secondary": "#52514e",
              "bar": "#2a78d6", "border": "black"},
    "dark": {"surface": "#1a1a19", "primary": "#ffffff", "secondary": "#c3c2b7",
             "bar": "#3987e5", "border": "#ffffff"},
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
    """Read the numeric columns of the `Rsquare.java` BED into an (n,5) array.
    An empty BED is legitimate (e.g. no call of a chromosome falls in a TR) and
    gives an empty array.
    """
    return np.loadtxt(data_lines(path), delimiter="\t", usecols=USECOLS,
                      ndmin=2).reshape(-1, len(USECOLS))


def select(data, code, args):
    """Keep the records of type `code` (of any type if `code` is None) that can
    be placed on both heatmaps, and count what was dropped and why."""
    keep = (np.ones(len(data), dtype=bool) if code is None
            else data[:, COL_TYPE].astype(np.int64) == code)
    n_of_type = int(keep.sum())
    if args.min_length is not None:
        keep &= data[:, COL_LENGTH] >= args.min_length
    if args.max_length is not None:
        keep &= data[:, COL_LENGTH] <= args.max_length
    n_filtered_out = n_of_type - int(keep.sum())

    uncomparable = keep & (data[:, COL_R2] < 0)  # R^2 == -1 sentinel
    keep &= ~uncomparable
    keep &= data[:, COL_AC] >= 1
    # A length of 0 cannot be placed on a log axis. `Rsquare.java` never emits
    # one (getLength is >=1 for every type it assigns), so this only guards a
    # hand-edited BED, but dropping the record from *every* panel is what keeps
    # the marginals the exact marginals of their heatmaps.
    n_nonpositive_length = int((keep & (data[:, COL_LENGTH] < 1)).sum())
    keep &= data[:, COL_LENGTH] >= 1
    return {
        "ac": data[keep, COL_AC].astype(np.int64),
        "length": data[keep, COL_LENGTH].astype(np.int64),
        "r2": np.clip(data[keep, COL_R2], 0.0, 1.0),
        "n_filtered_out": n_filtered_out,
        "n_uncomparable": int(uncomparable.sum()),
        "n_nonpositive_length": n_nonpositive_length,
    }


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


def style_spines(ax, theme):
    """Frame the panel with a border on all four sides."""
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(theme["border"])
        spine.set_linewidth(0.8)


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


def dump_counts(path, counts, x_edges, y_edges, x_lo_name, x_hi_name, extra=None):
    """Table view of exactly what one heatmap encodes: one row per x bin, then a
    final `all` row holding the column totals (the marginal of that stratum).
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
    ap.add_argument("inside_bed", help="output BED of Rsquare.java restricted to "
                                       "the calls inside TRs: "
                                       "CHROM,P,P+1,ID,type,length,AC,N,R^2,ACs,Ns,IDs")
    ap.add_argument("outside_bed", help="same, restricted to the calls outside TRs")
    ap.add_argument("total_samples", type=int,
                    help="T, total number of samples in the cohort; the x axis "
                         "of the top heatmaps is AF = AC/(2T)")
    ap.add_argument("svtype", type=str.upper, choices=SVTYPES,
                    help="plot only the records of this type, or of every "
                         "type with ALL (case-insensitive)")
    ap.add_argument("--out", help="output PNG (default: rsquare.<svtype>.png)")
    ap.add_argument("--x-bins", type=int, default=32,
                    help="target number of log-spaced bins on the x axis of "
                         "each heatmap, i.e. AC bins on the top row and length "
                         "bins on the bottom row (default 32)")
    ap.add_argument("--y-bins", type=int, default=20,
                    help="number of equal R^2 bins, in all panels (default 20)")
    ap.add_argument("--max-ac", type=int, help="upper AC limit (default: max in the data)")
    ap.add_argument("--min-length", type=float, help="keep only variants of length >= this")
    ap.add_argument("--max-length", type=float, help="keep only variants of length <= this")
    ap.add_argument("--color-scale", choices=("auto", "log", "linear"), default="auto",
                    help="count scale for the color ramp shared by every "
                         "heatmap (default auto: log when counts span >2 "
                         "orders of magnitude)")
    ap.add_argument("--no-median", action="store_true",
                    help="omit the per-bin median R^2 trend line of the "
                         "heatmaps (and the overall median in the marginals)")
    ap.add_argument("--dark", action="store_true", help="render on the dark surface")
    ap.add_argument("--dump-counts", action="store_true",
                    help="also write the bin counts of every heatmap as TSV "
                         "tables next to the PNG")
    ap.add_argument("--title-prefix", default="",
                    help="string prepended, followed by a space, to the title "
                         "of every column of panels (default: none)")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()
    if args.total_samples < 1:
        sys.exit("error: total_samples must be >= 1")
    n_alleles = 2 * args.total_samples  # AF = AC / (2T)
    code = None if args.svtype == SVTYPE_ALL else TYPE_CODES[args.svtype]

    # One column of panels per entry, left to right.
    strata = []
    for label, tag, path in (("TR regions", "inside_tr", args.inside_bed),
                             ("non-TR regions", "outside_tr", args.outside_bed)):
        data = load(path)
        # A wrong T (or a non-diploid callset) shows up here rather than as a
        # silently mislabelled axis.
        if len(data) and data[:, COL_AC].max() > n_alleles:
            sys.exit(f"error: {path} has AC up to {int(data[:, COL_AC].max()):,}, "
                     f"which exceeds 2T = {n_alleles:,} for T = "
                     f"{args.total_samples:,} samples")
        strata.append(dict(select(data, code, args), label=label, tag=tag))

    all_ac = np.concatenate([s["ac"] for s in strata])
    if all_ac.size == 0:
        sys.exit(f"error: no {args.svtype} records left to plot after filtering")
    max_ac = args.max_ac if args.max_ac else int(all_ac.max())
    n_over_max = 0
    for s in strata:
        in_range = s["ac"] <= max_ac
        n_over_max += int((~in_range).sum())
        for key in ("ac", "length", "r2"):
            s[key] = s[key][in_range]
    all_ac = np.concatenate([s["ac"] for s in strata])
    all_length = np.concatenate([s["length"] for s in strata])
    if all_ac.size == 0:
        sys.exit(f"error: no {args.svtype} records with AC <= {max_ac:,}")
    n_plotted = all_ac.size

    # ---- bin ---------------------------------------------------------------
    # Bin from the smallest value actually present, not from 1: Rsquare.java's
    # --min-af cutoff leaves nothing below ceil(2T*min_af), and its
    # --min-long-length cutoff leaves nothing below that length; starting at 1
    # would spend most of a log axis on empty bins.
    #
    # The edges come from both strata together and are then reused for each,
    # so the columns line up into one comparable grid: where a stratum has no
    # calls, its panel shows surface rather than a rebinned axis of its own.
    min_ac = max(1, int(all_ac.min()))
    x_edges = log_integer_edges(min_ac, max_ac, args.x_bins)
    y_edges = np.linspace(0.0, 1.0, args.y_bins + 1)
    min_length, max_length = max(1, int(all_length.min())), int(all_length.max())
    len_edges = log_integer_edges(min_length, max_length, args.x_bins)

    for s in strata:
        s["counts"] = np.histogram2d(s["ac"], s["r2"], bins=[x_edges, y_edges])[0]
        s["len_counts"] = np.histogram2d(s["length"], s["r2"],
                                         bins=[len_edges, y_edges])[0]

    theme = THEMES["dark" if args.dark else "light"]
    # matplotlib's default sequential colormap. .copy() so set_bad() does not
    # mutate the globally registered instance.
    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad(theme["surface"])  # empty cells take the surface color

    # One norm for every heatmap, so that equally dark cells mean equally many
    # calls anywhere in the figure.
    vmax = max(max(s["counts"].max(), s["len_counts"].max()) for s in strata)
    use_log = args.color_scale == "log" or (args.color_scale == "auto" and vmax >= 100)
    norm = LogNorm(vmin=1, vmax=max(vmax, 2)) if use_log else Normalize(vmin=0, vmax=vmax)
    # The marginals show each R^2 bin's fraction of its column's calls (each
    # marginal sums to 1). Each is scaled to its own largest fraction, since
    # they are read for their shape, not compared bar by bar.
    for s in strata:
        totals = s["counts"].sum(axis=0)
        s["r2_fractions"] = totals / totals.sum() if totals.sum() else totals

    # ---- draw --------------------------------------------------------------
    plt.rcParams.update({
        "font.size": 10,
        "text.color": theme["primary"],
        "axes.labelcolor": theme["secondary"],
        "xtick.color": theme["secondary"],
        "ytick.color": theme["secondary"],
    })
    column_width = PANEL + MARGINAL_GAP + MARGINAL
    col_x = (MARGIN_LEFT, MARGIN_LEFT + column_width + GUTTER)
    fig_w = col_x[1] + column_width + MARGIN_RIGHT
    fig_h = MARGIN_BOTTOM + 2 * PANEL + ROW_GAP + MARGIN_TOP
    row_y = (MARGIN_BOTTOM + PANEL + ROW_GAP, MARGIN_BOTTOM)  # AF row, length row
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=theme["surface"])

    def add_axes(x, y, w, h, **kwargs):
        ax = fig.add_axes((x / fig_w, y / fig_h, w / fig_w, h / fig_h), **kwargs)
        ax.set_facecolor(theme["surface"])
        return ax

    # Every panel shares the R^2 axis and every row shares its x axis, so a
    # band, a column or a cell lines up everywhere it appears.
    for c, stratum in enumerate(strata):
        first = strata[0].get("axes", {})
        ax = add_axes(col_x[c], row_y[0], PANEL, PANEL,
                      sharex=first.get("af"), sharey=first.get("af"))
        stratum["axes"] = {
            "af": ax,
            "dist": add_axes(col_x[c] + PANEL + MARGINAL_GAP, row_y[0], MARGINAL,
                             PANEL, sharey=ax),  # own x range: see r2_fractions
            "len": add_axes(col_x[c], row_y[1], PANEL, PANEL,
                            sharex=first.get("len"), sharey=ax),
        }

    # AC bins are built on integers (above); rescale to AF only for drawing, so
    # integer AC k occupies [(k-0.5)/2T, (k+0.5)/2T). Length is already in the
    # units of its axis and only needs the same half-unit centering.
    x_draw = (x_edges - 0.5) / n_alleles
    len_draw = len_edges - 0.5

    for c, stratum in enumerate(strata):
        ax, ax_dist, ax_len = (stratum["axes"][k] for k in ("af", "dist", "len"))
        left = c == 0
        n_calls = len(stratum["ac"])

        draw_heatmap(ax, x_draw, y_edges, stratum["counts"], cmap, norm, theme)
        draw_heatmap(ax_len, len_draw, y_edges, stratum["len_counts"], cmap, norm, theme)
        if n_calls == 0:
            # An empty stratum is a legitimate result (e.g. no INS inside TRs),
            # and an unexplained blank panel is not.
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
                    panel.plot(centers, medians, color=MEDIAN_COLOR, linewidth=1.0,
                               linestyle="-", zorder=4)

        # ---- top panel: R^2 versus AF --------------------------------------
        style_log_x(ax, x_draw, af_tick)
        ax.set_ylim(0.0, 1.0)
        ax.set_yticks(np.linspace(0, 1, 6))
        style_spines(ax, theme)
        ax.tick_params(length=3, width=0.8)
        ax.set_xlabel(f"AF")  #ax.set_xlabel(f"AF ({args.total_samples:,} samples)")

        # ---- marginal: the distribution of this column's R^2 values ---------
        # Summing the heatmap's columns rather than re-histogramming r2 keeps
        # the panels of a column describing the exact same calls, by
        # construction.
        ax_dist.barh(y_edges[:-1], stratum["r2_fractions"], height=np.diff(y_edges), align="edge",
                     color=theme["bar"], edgecolor=theme["surface"], linewidth=0.8,
                     zorder=2)
        median_all = float(np.median(stratum["r2"])) if n_calls else None
        if not args.no_median and n_calls:
            # The overall median, i.e. where this distribution's mass splits;
            # the heatmaps' lines are the same statistic taken per AF / per
            # length bin. Its value is printed in the column header, since the
            # marginal is too narrow to hold it.
            ax_dist.axhline(median_all, color=MEDIAN_COLOR, linewidth=1.0,
                            linestyle="-", zorder=4)
        ax_dist.set_xlim(0, (stratum["r2_fractions"].max() or 1) * 1.02)
        ax_dist.xaxis.set_visible(False)
        ax_dist.axis("off")  # no border: the bars are read only for their shape

        # ---- bottom panel: the same calls over call length -----------------
        style_log_x(ax_len, len_draw, count_tick)
        style_spines(ax_len, theme)
        ax_len.tick_params(length=3, width=0.8)
        ax_len.set_xlabel("Length (bp)")
        if left:
            ax.set_ylabel(Y_LABEL)
            ax_len.set_ylabel(Y_LABEL)

        # ---- header ---------------------------------------------------------
        header = " ".join(filter(None, (args.title_prefix, stratum["label"])))
        if not args.no_median and n_calls:
            header += f", median $R^2$ {median_all:.2f}"
        # Centered over the heatmap alone, not over heatmap and marginal.
        ax.text(0.5, HEADER_Y, header, transform=ax.transAxes, ha="center",
                va="bottom", color=theme["primary"], fontsize=10)

    out = args.out or f"rsquare.{args.svtype}.png"
    # Crop the figure to its outermost artists, with no surface border around
    # them. The pad is just wider than half a spine, which would otherwise be
    # cut off along the edge.
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", pad_inches=0.01,
                facecolor=theme["surface"])
    plt.close(fig)
    print(f"wrote {out}  ({n_plotted:,} {args.svtype}, {len(x_edges)-1} AF bins and "
          f"{len(len_edges)-1} length bins x {args.y_bins} R^2 bins, "
          f"{'log' if use_log else 'linear'} color scale up to {int(vmax):,} calls per bin)")
    print("  " + " · ".join(f"{s['label']}: {len(s['ac']):,}" for s in strata))
    # The figure has no room for what was dropped, so it is reported here.
    n_filtered_out = sum(s["n_filtered_out"] for s in strata)
    n_uncomparable = sum(s["n_uncomparable"] for s in strata)
    n_nonpositive_length = sum(s["n_nonpositive_length"] for s in strata)
    if n_filtered_out:
        print(f"  {n_filtered_out:,} excluded by filters")
    if n_uncomparable:
        print(f"  {n_uncomparable:,} with R^2=-1 (no comparable short) dropped")
    if n_nonpositive_length:
        print(f"  {n_nonpositive_length:,} with length <1 dropped")
    if n_over_max:
        print(f"  {n_over_max:,} above AF {max_ac / n_alleles:.3g} (AC {max_ac:,}) dropped")
    # Printing the ranges makes a T that does not match the BED obvious: with
    # Rsquare.java's --min-af cutoff, min AF should land on it.
    print(f"  AC {min_ac:,}..{max_ac:,}  ->  AF {min_ac/n_alleles:.5g}..{max_ac/n_alleles:.5g}"
          f"  (T = {args.total_samples:,}, 2T = {n_alleles:,})")
    print(f"  length {min_length:,}..{max_length:,} bp")

    if args.dump_counts:
        # Two tables per stratum, one per heatmap; each carries a final `all`
        # row, which is that stratum's marginal.
        stem = os.path.splitext(out)[0]
        for stratum in strata:
            tag = stratum["tag"]
            dump_counts(f"{stem}.{tag}.counts.tsv", stratum["counts"], x_edges,
                        y_edges, "ac_lo", "ac_hi",
                        extra={"af_lo": lambda k: f"{x_edges[0 if k is None else k]/n_alleles:.6g}",
                               "af_hi": lambda k: f"{(x_edges[-1 if k is None else k+1]-1)/n_alleles:.6g}"})
            dump_counts(f"{stem}.{tag}.length_counts.tsv", stratum["len_counts"],
                        len_edges, y_edges, "length_lo", "length_hi")


if __name__ == "__main__":
    main()
