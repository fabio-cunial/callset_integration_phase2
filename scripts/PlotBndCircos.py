#!/usr/bin/env python3

import argparse
import gzip
import math
import re
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path


MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
CONTIG_RE = re.compile(r"##contig=<ID=([^,>]+),length=(\d+)")


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


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


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


def read_vcf(path, nsamples_tag, min_nsamples):
    contig_lengths = {}
    records = []

    with open_text(path) as f:
        for line in f:
            if line.startswith("##contig="):
                m = CONTIG_RE.match(line.strip())
                if m:
                    contig_lengths[m.group(1)] = int(m.group(2))
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            alt = fields[4]
            info = parse_info_field(fields[7])

            if info.get("SVTYPE") and info.get("SVTYPE") != "BND":
                continue

            chrom2, pos2 = None, None
            if "CHR2" in info and "END" in info:
                chrom2 = info["CHR2"]
                try:
                    pos2 = int(info["END"])
                except ValueError:
                    pos2 = None
            if chrom2 is None or pos2 is None:
                chrom2, pos2 = parse_mate_from_alt(alt)

            if chrom2 is None or pos2 is None:
                continue

            # Keep only records whose two breakends map to standard chromosomes.
            if not is_standard_chrom(chrom) or not is_standard_chrom(chrom2):
                continue

            if nsamples_tag in info:
                try:
                    n_samples = float(info[nsamples_tag])
                except ValueError:
                    n_samples = 1.0
            else:
                n_samples = 1.0

            if min_nsamples is not None and n_samples < min_nsamples:
                continue

            records.append((chrom, pos, chrom2, pos2, n_samples))

    used_chroms = set()
    for chrom, _, chrom2, _, _ in records:
        used_chroms.add(chrom)
        used_chroms.add(chrom2)

    if not contig_lengths:
        max_seen = defaultdict(int)
        for chrom, pos, chrom2, pos2, _ in records:
            max_seen[chrom] = max(max_seen[chrom], pos)
            max_seen[chrom2] = max(max_seen[chrom2], pos2)
        contig_lengths = {k: max(v, 1) for k, v in max_seen.items()}
    else:
        # Keep only standard chromosomes that are actually referenced by retained records.
        contig_lengths = {
            k: v
            for k, v in contig_lengths.items()
            if k in used_chroms and is_standard_chrom(k)
        }

    return contig_lengths, records


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


def draw_plot(contig_lengths, records, output, nsamples_tag, title):
    chroms, layout = build_chrom_layout(contig_lengths, gap_rad=0.02)

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

    sample_values = [r[4] for r in records]
    if sample_values:
        smin = min(sample_values)
        smax = max(sample_values)
    else:
        smin, smax = 1.0, 1.0

    for chrom, pos, chrom2, pos2, n_samples in records:
        if chrom not in layout or chrom2 not in layout:
            continue
        x1, y1 = chrom_pos_to_xy(chrom, pos, contig_lengths, layout, outer_r - ring_w)
        x2, y2 = chrom_pos_to_xy(chrom2, pos2, contig_lengths, layout, outer_r - ring_w)

        if smax > smin:
            t = (n_samples - smin) / (smax - smin)
        else:
            t = 0.5
        lw = 0.3 + 3.2 * t
        alpha = 0.08 + 0.25 * t
        patch = curved_edge_patch(x1, y1, x2, y2, line_width=lw, alpha=alpha, color="#D1495B")
        ax.add_patch(patch)

    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.set_title(title, fontsize=13, pad=18)

    fig.text(0.02, 0.02, f"Edge width ~ INFO/{nsamples_tag}", fontsize=9)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Draw a Circos-like BND edge plot from a BND-only VCF(.gz)."
    )
    parser.add_argument("vcf", help="Input VCF or VCF.GZ containing BND records")
    parser.add_argument("output", help="Output image file (e.g. plot.png)")
    parser.add_argument(
        "--nsamples-tag",
        default="N_DISCOVERY_SAMPLES",
        help="INFO tag used to scale edge width (default: N_DISCOVERY_SAMPLES)",
    )
    parser.add_argument(
        "--title",
        default="BND Circos Plot",
        help="Plot title",
    )
    parser.add_argument(
        "--min-n-discovery-samples",
        type=float,
        default=None,
        help=(
            "Only plot records with INFO/N_DISCOVERY_SAMPLES >= this value "
            "(or >= selected --nsamples-tag value)"
        ),
    )
    args = parser.parse_args()

    contig_lengths, records = read_vcf(
        args.vcf,
        args.nsamples_tag,
        args.min_n_discovery_samples,
    )
    if not records:
        raise RuntimeError("No plottable BND records found in input VCF")

    draw_plot(contig_lengths, records, args.output, args.nsamples_tag, args.title)


if __name__ == "__main__":
    main()