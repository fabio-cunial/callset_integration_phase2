#!/usr/bin/env python3
# Warning: this script is AI-generated.
#
# Draws one Circos-like BND plot per sample, for the top-N samples of a
# multi-sample VCF. For every sample only the BNDs that are ALT in that sample
# are kept; then two BNDs of the same sample are collapsed iff both of their
# breakpoints are closer than a given distance (orientation is ignored, i.e.
# a BND and its mate, or two BNDs pointing in opposite directions, collapse).

import argparse
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
from array import array
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import PathPatch, Wedge  # noqa: E402
from matplotlib.path import Path  # noqa: E402


MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
CONTIG_RE = re.compile(r"##contig=<ID=([^,>]+),length=(\d+)")
INFO_ID_RE = re.compile(r"##INFO=<ID=([^,>]+)")
GT_SEP_RE = re.compile(r"[/|]")


def is_standard_chrom(chrom):
    if chrom is None:
        return False
    if not chrom.startswith("chr"):
        return False
    suffix = chrom[3:]
    if suffix.isdigit():
        n = int(suffix)
        return 1 <= n <= 22
    return suffix in {"X", "Y", "M", "MT"}


def natural_chrom_key(chrom):
    c = chrom
    if c.startswith("chr"):
        c = c[3:]
    if c.isdigit():
        return (0, int(c))
    if c == "X":
        return (1, 23)
    if c == "Y":
        return (1, 24)
    if c in ("M", "MT"):
        return (1, 25)
    return (2, c)


def open_text(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def have_bcftools():
    return shutil.which("bcftools") is not None


def parse_info_field(info):
    out = {}
    for token in info.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            out[k] = v
        elif token:
            out[token] = True
    return out


def parse_mate_from_alt(alt):
    m = MATE_RE.search(alt)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------


def read_header_lines(path, use_bcftools):
    if use_bcftools:
        out = subprocess.run(
            ["bcftools", "view", "-h", path],
            stdout=subprocess.PIPE,
            check=True,
            text=True,
        )
        return out.stdout.splitlines()
    lines = []
    with open_text(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            lines.append(line.rstrip("\n"))
    return lines


def parse_header(header_lines):
    """Returns (contig_lengths, samples, info_tags)."""
    contig_lengths = {}
    samples = []
    info_tags = set()
    for line in header_lines:
        if line.startswith("##contig="):
            m = CONTIG_RE.match(line)
            if m:
                contig_lengths[m.group(1)] = int(m.group(2))
        elif line.startswith("##INFO="):
            m = INFO_ID_RE.match(line)
            if m:
                info_tags.add(m.group(1))
        elif line.startswith("#CHROM"):
            fields = line.split("\t")
            if len(fields) > 9:
                samples = fields[9:]
    return contig_lengths, samples, info_tags


# ---------------------------------------------------------------------------
# Body iteration: both implementations yield
# (chrom, pos, alt_string, extra_tags_dict, genotype_strings)
# ---------------------------------------------------------------------------

BODY_TAGS = ("SVTYPE", "CHR2", "END")


def iter_body_bcftools(path, info_tags, nsamples_tag):
    wanted = [t for t in BODY_TAGS if t in info_tags]
    if nsamples_tag in info_tags and nsamples_tag not in wanted:
        wanted.append(nsamples_tag)

    columns = ["%CHROM", "%POS", "%ALT"] + ["%INFO/" + t for t in wanted]
    fmt = "\t".join(columns) + "[\t%GT]\n"
    n_fixed = len(columns)

    proc = subprocess.Popen(
        ["bcftools", "query", "-f", fmt, path],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1 << 20,
    )
    try:
        for line in proc.stdout:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < n_fixed:
                continue
            extra = dict(zip(wanted, fields[3:n_fixed]))
            yield fields[0], int(fields[1]), fields[2], extra, fields[n_fixed:]
    finally:
        proc.stdout.close()
        if proc.wait() != 0:
            raise RuntimeError("bcftools query failed on %s" % path)


def iter_body_python(path, info_tags, nsamples_tag):
    wanted = [t for t in BODY_TAGS] + [nsamples_tag]
    with open_text(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            info = parse_info_field(fields[7])
            extra = {}
            for tag in wanted:
                value = info.get(tag)
                if value is not None and value is not True:
                    extra[tag] = value
            # The VCF spec mandates GT as the first FORMAT subfield when present.
            gts = [s.split(":", 1)[0] for s in fields[9:]]
            yield fields[0], int(fields[1]), fields[4], extra, gts


def alt_alleles_of_genotype(gt):
    """Returns the set of 1-based ALT allele indexes carried by a genotype."""
    if not gt or gt == "." or gt == "./." or gt == ".|.":
        return ()
    out = set()
    for token in GT_SEP_RE.split(gt):
        if token.isdigit():
            k = int(token)
            if k > 0:
                out.add(k)
    return out


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def canonical_record(chrom1, pos1, chrom2, pos2):
    """Orientation-free representation: the two breakpoints in a fixed order."""
    if (natural_chrom_key(chrom1), pos1) <= (natural_chrom_key(chrom2), pos2):
        return (chrom1, pos1, chrom2, pos2)
    return (chrom2, pos2, chrom1, pos1)


def load_bnds(path, use_bcftools, info_tags, samples, nsamples_tag, min_ns, max_ns):
    """Returns (records, per_sample), where records is a list of canonical BNDs
    and per_sample maps a sample name to the array of its ALT record indexes."""
    records = []
    per_sample = {s: array("i") for s in samples}
    iterator = iter_body_bcftools if use_bcftools else iter_body_python

    n_lines = 0
    for chrom, pos, alt_string, extra, gts in iterator(path, info_tags, nsamples_tag):
        n_lines += 1
        svtype = extra.get("SVTYPE")
        if svtype and svtype != "." and svtype != "BND":
            continue
        if not is_standard_chrom(chrom):
            continue

        value = extra.get(nsamples_tag)
        if value is not None and value != ".":
            try:
                n_discovery = float(value)
            except ValueError:
                n_discovery = 1.0
            if min_ns is not None and n_discovery < min_ns:
                continue
            if max_ns is not None and n_discovery > max_ns:
                continue

        alts = alt_string.split(",")
        # The CHR2/END fallback is only unambiguous for bi-allelic records.
        chr2_info, end_info = None, None
        if len(alts) == 1 and extra.get("CHR2", ".") != "." and extra.get("END", ".") != ".":
            chr2_info = extra["CHR2"]
            try:
                end_info = int(extra["END"])
            except ValueError:
                end_info = None

        # One entry of `records` per ALT allele that is a usable BND.
        allele_to_index = {}
        for k, alt in enumerate(alts, start=1):
            chrom2, pos2 = parse_mate_from_alt(alt)
            if chrom2 is None and chr2_info is not None and end_info is not None:
                chrom2, pos2 = chr2_info, end_info
            if chrom2 is None or pos2 is None:
                continue
            if not is_standard_chrom(chrom2):
                continue
            allele_to_index[k] = len(records)
            records.append(canonical_record(chrom, pos, chrom2, pos2))
        if not allele_to_index:
            continue

        for i, gt in enumerate(gts):
            if i >= len(samples):
                break
            alleles = alt_alleles_of_genotype(gt)
            if not alleles:
                continue
            dest = per_sample[samples[i]]
            for k in alleles:
                index = allele_to_index.get(k)
                if index is not None:
                    dest.append(index)

    sys.stderr.write("Scanned %d VCF records; kept %d BND alleles\n" % (n_lines, len(records)))
    return records, per_sample


# ---------------------------------------------------------------------------
# Collapsing
# ---------------------------------------------------------------------------


def collapse_records(indexes, records, distance):
    """Single-linkage collapse of the BNDs of one sample: two BNDs are merged
    iff both of their breakpoints are within `distance` of each other. Returns
    a list of (chrom1, pos1, chrom2, pos2, n_collapsed)."""
    n = len(indexes)
    if n == 0:
        return []

    parent = list(range(n))

    def find(x):
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    # Bin the breakpoint pairs on a `distance`-wide grid: two BNDs that must be
    # merged always fall into adjacent bins, so only 9 bins have to be probed.
    bins = defaultdict(list)
    for i, index in enumerate(indexes):
        c1, p1, c2, p2 = records[index]
        bins[(c1, c2, p1 // distance, p2 // distance)].append(i)

    for i, index in enumerate(indexes):
        c1, p1, c2, p2 = records[index]
        b1, b2 = p1 // distance, p2 // distance
        for d1 in (-1, 0, 1):
            for d2 in (-1, 0, 1):
                for j in bins.get((c1, c2, b1 + d1, b2 + d2), ()):
                    if j <= i:
                        continue
                    _, q1, _, q2 = records[indexes[j]]
                    if abs(p1 - q1) <= distance and abs(p2 - q2) <= distance:
                        union(i, j)

    clusters = defaultdict(list)
    for i in range(n):
        clusters[find(i)].append(i)

    out = []
    for members in clusters.values():
        c1, _, c2, _ = records[indexes[members[0]]]
        p1 = sum(records[indexes[m]][1] for m in members) // len(members)
        p2 = sum(records[indexes[m]][3] for m in members) // len(members)
        out.append((c1, p1, c2, p2, float(len(members))))
    return out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def build_chrom_layout(contig_lengths, gap_rad):
    chroms = sorted(contig_lengths.keys(), key=natural_chrom_key)
    total_bp = sum(contig_lengths[c] for c in chroms)
    if total_bp <= 0:
        raise ValueError("No contig length information available")

    total_gap = gap_rad * len(chroms)
    usable = 2.0 * math.pi - total_gap
    if usable <= 0:
        raise ValueError("Too many chromosomes for chosen gap")

    layout = {}
    angle = 0.0
    for chrom in chroms:
        span = usable * (contig_lengths[chrom] / total_bp)
        layout[chrom] = (angle, angle + span)
        angle += span + gap_rad
    return chroms, layout


def chrom_pos_to_xy(chrom, pos, contig_lengths, layout, radius):
    start, end = layout[chrom]
    clen = max(contig_lengths[chrom], 1)
    t = min(max((pos - 1) / clen, 0.0), 1.0)
    theta = start + t * (end - start)
    x = radius * math.cos(theta)
    y = radius * math.sin(theta)
    return x, y


def curved_edge_patch(x1, y1, x2, y2, line_width, alpha, color):
    ctrl_scale = 0.15
    cx1 = x1 * ctrl_scale
    cy1 = y1 * ctrl_scale
    cx2 = x2 * ctrl_scale
    cy2 = y2 * ctrl_scale

    path = Path(
        [(x1, y1), (cx1, cy1), (cx2, cy2), (x2, y2)],
        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4],
    )
    return PathPatch(path, lw=line_width, edgecolor=color, facecolor="none", alpha=alpha)


def draw_plot(contig_lengths, chroms, layout, records, output, title, footer):
    fig = plt.figure(figsize=(12, 12), dpi=180)
    ax = fig.add_subplot(111)
    ax.set_aspect("equal")
    ax.axis("off")

    outer_r = 1.0
    ring_w = 0.06

    palette = [
        "#2D6A4F",
        "#40916C",
        "#1B4332",
        "#A7C957",
        "#52796F",
        "#4D908E",
        "#76C893",
        "#34A0A4",
    ]

    for i, chrom in enumerate(chroms):
        start, end = layout[chrom]
        theta1 = math.degrees(start)
        theta2 = math.degrees(end)
        w = Wedge(
            (0, 0),
            r=outer_r,
            theta1=theta1,
            theta2=theta2,
            width=ring_w,
            facecolor=palette[i % len(palette)],
            edgecolor="white",
            lw=0.6,
            alpha=0.95,
        )
        ax.add_patch(w)

        mid = 0.5 * (start + end)
        tx = (outer_r + 0.08) * math.cos(mid)
        ty = (outer_r + 0.08) * math.sin(mid)
        rot = math.degrees(mid)
        label = chrom[3:] if chrom.startswith("chr") else chrom
        ax.text(
            tx,
            ty,
            label,
            ha="center",
            va="center",
            fontsize=8,
            rotation=rot if -90 <= rot <= 90 else rot + 180,
            rotation_mode="anchor",
        )

    weights = [r[4] for r in records]
    if weights:
        smin = min(weights)
        smax = max(weights)
    else:
        smin, smax = 1.0, 1.0

    for chrom, pos, chrom2, pos2, weight in records:
        if chrom not in layout or chrom2 not in layout:
            continue
        x1, y1 = chrom_pos_to_xy(chrom, pos, contig_lengths, layout, outer_r - ring_w)
        x2, y2 = chrom_pos_to_xy(chrom2, pos2, contig_lengths, layout, outer_r - ring_w)

        if smax > smin:
            t = (weight - smin) / (smax - smin)
        else:
            t = 0.5
        lw = 0.3 + 3.2 * t
        alpha = 0.08 + 0.25 * t
        patch = curved_edge_patch(x1, y1, x2, y2, line_width=lw, alpha=alpha, color="#D1495B")
        ax.add_patch(patch)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_title(title, fontsize=13, pad=18)

    fig.text(0.02, 0.02, footer, fontsize=9)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------


def sanitize(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Draw one Circos-like BND plot per sample, for the top-N samples of a "
            "multi-sample VCF(.gz)."
        )
    )
    parser.add_argument("vcf", help="Input multi-sample VCF, VCF.GZ or BCF containing BND records")
    parser.add_argument("outdir", help="Output directory; one <sample>.png is written per sample")
    parser.add_argument(
        "--collapse-distance",
        type=int,
        default=1000,
        help=(
            "Two BNDs of the same sample are collapsed iff both of their breakpoints "
            "are within this many bp of each other, regardless of orientation "
            "(default: 1000)"
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of samples to plot, ranked by collapsed BND count (default: 10)",
    )
    parser.add_argument(
        "--nsamples-tag",
        default="N_DISCOVERY_SAMPLES",
        help=(
            "INFO tag used by --min-n-discovery-samples/--max-n-discovery-samples "
            "(default: N_DISCOVERY_SAMPLES)"
        ),
    )
    parser.add_argument(
        "--title",
        default="BND Circos Plot",
        help="Plot title prefix; the sample name is appended",
    )
    parser.add_argument(
        "--min-n-discovery-samples",
        type=float,
        default=None,
        help="Only consider records with INFO/<--nsamples-tag> >= this value",
    )
    parser.add_argument(
        "--max-n-discovery-samples",
        type=float,
        default=None,
        help="Only consider records with INFO/<--nsamples-tag> <= this value",
    )
    parser.add_argument(
        "--no-bcftools",
        action="store_true",
        help="Parse the VCF in pure Python instead of shelling out to bcftools",
    )
    args = parser.parse_args()

    if args.collapse_distance < 1:
        parser.error("--collapse-distance must be at least 1")

    use_bcftools = (not args.no_bcftools) and have_bcftools()
    if not use_bcftools and not args.no_bcftools:
        sys.stderr.write("bcftools not found in PATH: falling back to the Python parser\n")

    contig_lengths, samples, info_tags = parse_header(read_header_lines(args.vcf, use_bcftools))
    if not samples:
        raise RuntimeError("Input VCF has no samples")

    records, per_sample = load_bnds(
        args.vcf,
        use_bcftools,
        info_tags,
        samples,
        args.nsamples_tag,
        args.min_n_discovery_samples,
        args.max_n_discovery_samples,
    )
    if not records:
        raise RuntimeError("No plottable BND records found in input VCF")

    # Collapse every sample independently and rank samples by collapsed count.
    collapsed = {}
    for sample in samples:
        indexes = per_sample[sample]
        if not len(indexes):
            continue
        collapsed[sample] = collapse_records(indexes, records, args.collapse_distance)
    if not collapsed:
        raise RuntimeError("No sample carries any plottable BND")

    ranking = sorted(collapsed.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    top = ranking[: args.top_n]

    # A single layout, shared by all plots, so that they are comparable.
    if contig_lengths:
        layout_lengths = {k: v for k, v in contig_lengths.items() if is_standard_chrom(k)}
    else:
        layout_lengths = {}
    if not layout_lengths:
        max_seen = defaultdict(int)
        for c1, p1, c2, p2 in records:
            max_seen[c1] = max(max_seen[c1], p1)
            max_seen[c2] = max(max_seen[c2], p2)
        layout_lengths = {k: max(v, 1) for k, v in max_seen.items()}
    chroms, layout = build_chrom_layout(layout_lengths, gap_rad=0.02)

    title_filters = []
    if args.min_n_discovery_samples is not None:
        title_filters.append(f"min {args.nsamples_tag}>={args.min_n_discovery_samples:g}")
    if args.max_n_discovery_samples is not None:
        title_filters.append(f"max {args.nsamples_tag}<={args.max_n_discovery_samples:g}")

    os.makedirs(args.outdir, exist_ok=True)
    for sample, sample_records in top:
        n_raw = len(per_sample[sample])
        title = f"{args.title} - {sample}"
        if title_filters:
            title = f"{title} ({'; '.join(title_filters)})"
        footer = (
            f"{len(sample_records)} BNDs after collapsing {n_raw} ALT BNDs at "
            f"{args.collapse_distance} bp; edge width ~ number of collapsed BNDs"
        )
        output = os.path.join(args.outdir, sanitize(sample) + ".png")
        draw_plot(layout_lengths, chroms, layout, sample_records, output, title, footer)
        sys.stderr.write("%s: %d collapsed BNDs (%d raw) -> %s\n" % (sample, len(sample_records), n_raw, output))


if __name__ == "__main__":
    main()
