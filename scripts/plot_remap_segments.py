#!/usr/bin/env python3
"""
Plots the `truvari anno remap` annotations of every INS record in every VCF of a
directory, as one horizontal track per record, all of the same width.

For a record, the ALT sequence is drawn as a horizontal band of fixed width
(every record is normalized to the same width, regardless of its ALT length).
Every alignment in `remap_segments` is drawn as a thick arrow that spans the
alignment's interval in ALT, and that points right (+) or left (-) according to
the alignment's orientation. Arrows are filled with a gradient: a fixed rainbow
gradient is assigned to the whole `remap_coords` interval, and every alignment
gets the subgradient of its own reference interval, reversed if the orientation
is `-`. Portions of ALT that are not covered by any alignment are drawn as a
thin horizontal line.

Rows are clustered by architecture before being drawn. Since both axes are
normalized (ALT to [0,1], reference to a fraction of the `remap_coords`
interval), the similarity is locus- and size-invariant: records with the same
structural pattern cluster together even if they lie on different chromosomes
and differ in length by orders of magnitude. Two measures are available:

  --metric edit (default)
      The `remap_coords` interval is partitioned by `--ref-splits` splitpoints,
      and every pair of splitpoints is a character of an alphabet; every
      alignment becomes the character of the pair closest to its reference
      interval, primed if the alignment is reverse-complemented. ALT is
      partitioned by `--alt-chunks` splitpoints, and every quantized length is
      a character of a second, disjoint alphabet; the gaps between consecutive
      alignments become characters of this alphabet. A record is then the
      string of its alignments and gaps, in ALT order, and the distance is the
      edit distance between two strings, divided by the length of the longer
      one. See `record_tokens()` and `edit_distance_matrix()`.

      The quantization is what absorbs the variation of the breakpoints
      between two records of the same architecture: two alignments are the
      same character only if they fall in the same pair of splitpoints, and a
      pair of records whose breakpoints differ by more than a splitpoint pays
      a full replacement. `--ref-splits` should therefore stay coarser than
      that variation (20 splitpoints tolerate breakpoints that move by up to
      ~2.5% of the `remap_coords` interval). A finer partition can be used
      together with `--sub-cost graded`, which charges a replacement in
      proportion to how far apart the two characters are, and which is
      insensitive to the granularity.

  --metric profile
      Mean per-bin disagreement between the two rows, rasterized along ALT.
      See `record_profiles()` and `profile_distance_matrix()`.

Usage:
    python plot_remap_segments.py <input_dir> [-o out.png]
"""

import argparse
import glob
import gzip
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager

import matplotlib

matplotlib.use("Agg")  # Save to disk without displaying.

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon, Rectangle

try:
    from scipy.cluster.hierarchy import fcluster, leaves_list, linkage, optimal_leaf_ordering
    from scipy.spatial.distance import squareform

    HAS_SCIPY = True
except ImportError:  # Falls back to the pure-numpy implementations below.
    HAS_SCIPY = False

try:  # Computes the whole matrix of edit distances in C, on every core.
    from rapidfuzz import process as rapidfuzz_process
    from rapidfuzz.distance import Levenshtein as rapidfuzz_levenshtein

    HAS_RAPIDFUZZ = True
except ImportError:  # Falls back to the pure-python implementation below.
    HAS_RAPIDFUZZ = False

COORDS_RE = re.compile(r"^(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)$")
SEGMENT_RE = re.compile(
    r"^(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)\((?P<ori>[+-])\):q(?P<qstart>\d+)-(?P<qend>\d+)$"
)
SEQUENCE_RE = re.compile(r"^[ACGTNacgtn]+$")

# Arrows narrower than this fraction of the track width are widened, so that
# short alignments stay visible.
MIN_ARROW_WIDTH = 0.002
# Arrowhead length, as a fraction of the track width and of the arrow length.
MAX_HEAD_WIDTH = 0.012
MAX_HEAD_FRACTION = 0.35
# Number of orders of magnitude of ALT length over which the (small) length
# term of the distance saturates.
LENGTH_DECADES = 2.0
# Keeps the alphabet of the gaps disjoint from the alphabet of the alignments.
GAP_ALPHABET_BASE = 1 << 40
# Resolution, in decades, at which the ALT length enters the distance. Records
# that agree on their string and on their length bucket are indistinguishable,
# and are collapsed into one key.
LENGTH_QUANTUM = 0.02
# Rows of the matrix of the distances that are materialized at a time. Small
# enough that the per-thread blocks stay modest on the largest cohorts.
KEY_CHUNK = 256
# Axes computed by --embed pca, shown as two panels of two axes each.
PCA_DIMS = 4
# Where the tree of --source is cut into flat clusters, as a fraction of its
# tallest merge. The Euclidean distance between distributions has no absolute
# scale to compare against, so the cut has to be relative.
SOURCE_THRESHOLD_FRACTION = 0.3
# Default cuts of the tree into flat clusters, one per metric.
EDIT_THRESHOLD = 0.30
PROFILE_THRESHOLD = 0.12
# The two orientations, wherever they are drawn as two series. Checked for
# colorblind separation against each other and for contrast against white.
FORWARD_COLOR = "#2a78d6"
REVERSE_COLOR = "#eb6834"
# Range of the inset of the panel of the lengths of the alignments, in bp. The
# panel itself spans the whole range, where a handful of very long alignments
# compress everything else into the leftmost bins; the inset is that left part,
# re-binned over its own range.
INSET_MAX_LENGTH = 10000


TIMINGS = []


@contextmanager
def phase(name):
    """
    Times a stage of the script. The breakdown is printed at the end: several
    stages (parsing, scipy's linkage and leaf ordering) are single-threaded by
    nature, so knowing which one dominates is the only way to tell whether a
    slow run can be helped by more cores at all.
    """
    start = time.time()
    yield
    TIMINGS.append((name, time.time() - start))


# ------------------------------- Parsing ------------------------------------

def open_vcf(path):
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(info):
    out = {}
    for field in info.split(";"):
        if not field:
            continue
        key, _, value = field.partition("=")
        out[key] = value
    return out


def parse_coords(coords):
    match = COORDS_RE.match(coords)
    if match is None:
        raise ValueError("Invalid remap_coords value: %s" % coords)
    return match.group("chrom"), int(match.group("start")), int(match.group("end"))


def parse_segments(segments):
    """
    Returns a list of (chrom,ref_start,ref_end,orientation,alt_start,alt_end).
    """
    out = []
    for piece in segments.strip().strip("|").split("|"):
        if not piece:
            continue
        match = SEGMENT_RE.match(piece)
        if match is None:
            raise ValueError("Invalid remap_segments alignment: %s" % piece)
        out.append(
            (
                match.group("chrom"),
                int(match.group("start")),
                int(match.group("end")),
                match.group("ori"),
                int(match.group("qstart")),
                int(match.group("qend")),
            )
        )
    return out


def load_records(path):
    """
    Returns a list of dictionaries, one per INS record with remap annotations.
    """
    out = []
    with open_vcf(path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            tokens = line.rstrip("\n").split("\t")
            if len(tokens) < 8:
                continue
            info = parse_info(tokens[7])
            if "remap_coords" not in info or "remap_segments" not in info:
                continue
            try:
                coords = parse_coords(info["remap_coords"])
                segments = parse_segments(info["remap_segments"])
            except ValueError as error:
                print("Skipping %s:%s - %s" % (tokens[0], tokens[1], error), file=sys.stderr)
                continue
            if not segments:
                continue
            # The ALT sequence length, used to normalize every track to the same
            # width. Falls back to the rightmost aligned ALT position for
            # symbolic ALTs.
            alt = tokens[4]
            if SEQUENCE_RE.match(alt) is not None:
                alt_length = len(alt)
            else:
                alt_length = max(segment[5] for segment in segments)
            alt_length = max(alt_length, max(segment[5] for segment in segments))
            out.append(
                {
                    "file": os.path.basename(path),
                    "chrom": tokens[0],
                    "pos": tokens[1],
                    "id": tokens[2],
                    "coords": coords,
                    "ori": info.get("remap_ori", ""),
                    "segments": segments,
                    "alt_length": alt_length,
                }
            )
    return out


# ------------------------------ Similarity ----------------------------------

def color_fractions(record, seg_start, seg_end):
    """
    Position of a reference interval inside the `remap_coords` interval, as a
    pair of fractions in [0,1]. This is exactly the subgradient of the segment.
    """
    _, ref_start, ref_end = record["coords"]
    ref_span = max(ref_end - ref_start, 1)
    return (
        float(np.clip((seg_start - ref_start) / ref_span, 0, 1)),
        float(np.clip((seg_end - ref_start) / ref_span, 0, 1)),
    )


def record_profiles(records, n_bins):
    """
    Rasterizes every record into a fixed-length profile of the reference
    fraction (i.e. of the gradient value) and of the orientation, sampled at
    `n_bins` equally spaced points along ALT. This is what `--metric profile`
    is computed on: it is what the eye sees, with both axes normalized, so it
    does not depend on the locus or on the size of the event.

    Returns (colors,orientations,covered), each of shape (len(records),n_bins).
    Bins that no alignment covers are marked in `covered` and ignored by
    `profile_distance_matrix()`.
    """
    n = len(records)
    colors = np.zeros((n, n_bins))
    orientations = np.zeros((n, n_bins), dtype=np.int8)
    covered = np.zeros((n, n_bins), dtype=bool)
    centers = (np.arange(n_bins) + 0.5) / n_bins
    for i, record in enumerate(records):
        alt_length = max(record["alt_length"], 1)
        positions = centers * alt_length
        # Longest alignments are written last, so that they win wherever
        # several alignments overlap in ALT.
        for _, seg_start, seg_end, ori, alt_start, alt_end in sorted(
            record["segments"], key=lambda segment: segment[5] - segment[4]
        ):
            left, right = alt_start - 1, alt_end
            inside = (positions >= left) & (positions < right)
            if not inside.any():
                continue
            start_fraction, end_fraction = color_fractions(record, seg_start, seg_end)
            u = (positions[inside] - left) / max(right - left, 1)
            if ori == "+":
                colors[i, inside] = start_fraction + u * (end_fraction - start_fraction)
            else:  # Reversed subgradient, as in the drawing.
                colors[i, inside] = end_fraction - u * (end_fraction - start_fraction)
            orientations[i, inside] = 1 if ori == "+" else -1
            covered[i, inside] = True
    return colors, orientations, covered


def profile_distance_matrix(records, colors, orientations, covered, w_ori, gap, len_weight):
    """
    Mean per-bin disagreement between two rows, in [0,1]:

      - both bins uncovered:        0
      - exactly one bin covered:    `gap`
      - both covered:               (1-w_ori)*|dcolor| + w_ori*[different strand]
    """
    n, n_bins = colors.shape
    w_color = 1.0 - w_ori
    out = np.zeros((n, n))
    for i in range(n):
        both = covered[i] & covered
        one = covered[i] ^ covered
        color = np.abs(colors[i] - colors) * both
        strand = (orientations[i] != orientations) & both
        out[i] = (w_color * color.sum(1) + w_ori * strand.sum(1) + gap * one.sum(1)) / n_bins
    out = (out + out.T) / 2  # Guards against asymmetry from rounding.
    np.fill_diagonal(out, 0.0)
    return add_length_term(out, records, len_weight)


def add_length_term(distances, records, len_weight):
    """
    Blends into `distances` a small term on the ratio of the ALT lengths, so
    that size still breaks ties between records of identical architecture.
    Keeps the result in [0,1].
    """
    if len_weight <= 0:
        return distances
    lengths = np.log10(np.maximum([record["alt_length"] for record in records], 1))
    ratio = np.minimum(np.abs(lengths[:, None] - lengths[None, :]) / LENGTH_DECADES, 1.0)
    return (distances + len_weight * ratio) / (1.0 + len_weight)


def record_tokens(record, n_ref_splits, n_alt_chunks):
    """
    Encodes a record as a string over two disjoint alphabets:

      - an alignment is the character of the pair of splitpoints of
        `remap_coords` closest to its reference interval, primed (a distinct
        character) when the alignment is reverse-complemented;
      - a maximal region of ALT that no alignment covers is the character of
        its length, quantized to `n_alt_chunks` chunks of ALT.

    Characters are emitted in ALT order, and are returned as integers. Gaps
    shorter than half a chunk quantize to zero and are dropped, which is what
    makes the encoding robust to small differences in the alignments.
    """
    _, ref_start, ref_end = record["coords"]
    ref_span = max(ref_end - ref_start, 1)
    alt_length = max(record["alt_length"], 1)
    out = []

    def push_gap(length):
        chunks = int(round(length / alt_length * n_alt_chunks))
        if chunks >= 1:
            out.append(GAP_ALPHABET_BASE + chunks)

    cursor = 0
    for _, seg_start, seg_end, ori, alt_start, alt_end in sorted(
        record["segments"], key=lambda segment: (segment[4], segment[5])
    ):
        left, right = alt_start - 1, alt_end
        if left > cursor:
            push_gap(left - cursor)
        start_split = int(round((seg_start - ref_start) / ref_span * n_ref_splits))
        end_split = int(round((seg_end - ref_start) / ref_span * n_ref_splits))
        start_split = min(max(start_split, 0), n_ref_splits)
        end_split = min(max(end_split, 0), n_ref_splits)
        if start_split > end_split:
            start_split, end_split = end_split, start_split
        out.append(
            (start_split * (n_ref_splits + 1) + end_split) * 2 + (0 if ori == "+" else 1)
        )
        cursor = max(cursor, right)
    if cursor < alt_length:
        push_gap(alt_length - cursor)
    return tuple(out)


def decode_token(token, n_ref_splits):
    """
    Inverse of the encoding of `record_tokens()`. Returns
    (is_gap,first,second,strand): (chunks,0,0) for a gap character, and
    (start_split,end_split,strand) for an alignment character.
    """
    if token >= GAP_ALPHABET_BASE:
        return True, token - GAP_ALPHABET_BASE, 0, 0
    strand = token & 1
    pair = token >> 1
    return False, pair // (n_ref_splits + 1), pair % (n_ref_splits + 1), strand


def substitution_costs(alphabet, n_ref_splits, n_alt_chunks, graded, w_ori):
    """
    Cost of replacing one character by another, for every pair of characters
    that occurs in the data. With `graded`, replacing a character by a nearby
    one (a slightly different reference interval, or a slightly different gap
    length) costs less than replacing it by a distant one, which keeps records
    that quantize just across a splitpoint from looking unrelated. Otherwise
    every replacement of distinct characters costs 1, as in plain edit
    distance.
    """
    if not graded:
        return None
    decoded = {token: decode_token(token, n_ref_splits) for token in alphabet}
    out = {}
    for x, (x_gap, x_first, x_second, x_strand) in decoded.items():
        for y, (y_gap, y_first, y_second, y_strand) in decoded.items():
            if x == y:
                continue
            if x_gap != y_gap:  # A gap never resembles an alignment.
                out[(x, y)] = 1.0
            elif x_gap:
                out[(x, y)] = min(abs(x_first - y_first) / max(n_alt_chunks, 1), 1.0)
            else:
                shift = (abs(x_first - y_first) + abs(x_second - y_second)) / (
                    2.0 * max(n_ref_splits, 1)
                )
                out[(x, y)] = min(
                    (1.0 - w_ori) * min(shift, 1.0) + w_ori * (x_strand != y_strand), 1.0
                )
    return out


def edit_distance(x, y, costs):
    """
    Levenshtein distance between two strings of characters, with unit
    insertion and deletion costs, and with the replacement costs of `costs`
    (unit costs if `costs` is None).
    """
    if x == y:
        return 0.0
    if not x:
        return float(len(y))
    if not y:
        return float(len(x))
    previous = [float(j) for j in range(len(y) + 1)]
    current = [0.0] * (len(y) + 1)
    for i in range(1, len(x) + 1):
        current[0] = float(i)
        x_token = x[i - 1]
        for j in range(1, len(y) + 1):
            y_token = y[j - 1]
            if x_token == y_token:
                replace = previous[j - 1]
            elif costs is None:
                replace = previous[j - 1] + 1.0
            else:
                replace = previous[j - 1] + costs[(x_token, y_token)]
            current[j] = min(previous[j] + 1.0, current[j - 1] + 1.0, replace)
        previous, current = current, previous
    return previous[len(y)]


def all_pairs_edit_distances(keys, costs):
    """
    Matrix of the edit distances between every pair of the strings `keys`.
    Unit replacement costs go to rapidfuzz, which runs the bit-parallel
    algorithm in C on every core; graded replacement costs have no equivalent
    there, and fall back to the pure-python dynamic program, which is orders
    of magnitude slower on large inputs.
    """
    if costs is None and HAS_RAPIDFUZZ:
        return rapidfuzz_process.cdist(
            keys,
            keys,
            scorer=rapidfuzz_levenshtein.distance,
            dtype=np.float32,
            workers=-1,
        )
    out = np.zeros((len(keys), len(keys)), dtype=np.float32)
    for a in range(len(keys)):
        for b in range(a + 1, len(keys)):
            out[a, b] = out[b, a] = edit_distance(keys[a], keys[b], costs)
    return out


def record_key(record, args):
    """
    Everything the distance of a record depends on: its string of characters,
    and the bucket of its ALT length. Two records with the same key are at
    distance 0 from each other and equidistant from every other record, so the
    whole pipeline can work on the distinct keys instead of on the records,
    which is what makes large cohorts tractable.
    """
    tokens = record_tokens(record, args.ref_splits, args.alt_chunks)
    if args.len_weight <= 0:  # Length ignored altogether.
        return tokens, 0
    decades = np.log10(max(record["alt_length"], 1))
    return tokens, int(round(decades / LENGTH_QUANTUM))


def edit_distance_matrix(records, args):
    """
    Edit distance between the strings of `record_tokens()`, divided by the
    length of the longer string, so that the result lies in [0,1] and does not
    grow with the number of alignments, plus the usual term on the ALT length.

    Returns (distances,index,counts), where `distances` is over the distinct
    keys of `record_key()`, `index` maps a record to its key, and `counts` is
    the number of records per key.
    """
    unique, index = {}, []
    for record in records:
        index.append(unique.setdefault(record_key(record, args), len(unique)))
    keys = list(unique.keys())
    # Several keys can share the same string, when they differ only by length.
    strings, string_of_key = {}, []
    for tokens, _ in keys:
        string_of_key.append(strings.setdefault(tokens, len(strings)))
    string_of_key = np.array(string_of_key)
    distinct = list(strings.keys())

    alphabet = {token for string in distinct for token in string}
    costs = substitution_costs(
        alphabet, args.ref_splits, args.alt_chunks, args.sub_cost == "graded", args.w_orientation
    )
    print(
        "%d records, %d distinct keys, %d distinct strings over %d characters (%s)"
        % (
            len(records),
            len(keys),
            len(distinct),
            len(alphabet),
            "rapidfuzz" if costs is None and HAS_RAPIDFUZZ else "pure python",
        ),
        file=sys.stderr,
    )
    between_strings = all_pairs_edit_distances(distinct, costs)

    # Expands to the keys in chunks. The square matrix is the memory
    # bottleneck of the whole script, so the normalization by the longer
    # string and the term on the ALT length are applied per chunk: computed
    # over the full matrix, either would allocate a second copy of it.
    n_keys = len(keys)
    out = np.empty((n_keys, n_keys), dtype=np.float32)
    buckets = np.array([bucket for _, bucket in keys], dtype=np.float32) * LENGTH_QUANTUM
    lengths = np.array(
        [max(len(distinct[string]), 1) for string in string_of_key], dtype=np.float32
    )
    threads = thread_count(args.threads)

    def work(start, stop):
        block = between_strings[np.ix_(string_of_key[start:stop], string_of_key)]
        block /= np.maximum(lengths[start:stop, None], lengths[None, :])
        if args.len_weight > 0:
            ratio = np.minimum(
                np.abs(buckets[start:stop, None] - buckets[None, :]) / LENGTH_DECADES, 1.0
            )
            block = (block + args.len_weight * ratio) / (1.0 + args.len_weight)
        out[start:stop] = block

    parallel_chunks(n_keys, KEY_CHUNK, work, threads)
    if n_keys > 10000:  # The square matrix, plus what the clustering needs.
        print(
            "%d distinct keys need about %.1f GB: lower --ref-splits or --alt-chunks, or "
            "--max-records, if that is too much"
            % (n_keys, 3.0 * n_keys * n_keys * 4 / 1e9),
            file=sys.stderr,
        )
    return out, np.array(index), np.bincount(index, minlength=n_keys)


def average_linkage(distances):
    """
    Pure-numpy average linkage, in scipy's linkage-matrix format. Used only
    when scipy is unavailable.
    """
    n = distances.shape[0]
    current = distances.astype(float).copy()
    np.fill_diagonal(current, np.inf)
    sizes = np.ones(n)
    ids = np.arange(n)
    active = np.ones(n, dtype=bool)
    out = np.zeros((n - 1, 4))
    for m in range(n - 1):
        masked = np.where(active[:, None] & active[None, :], current, np.inf)
        i, j = np.unravel_index(np.argmin(masked), masked.shape)
        if i > j:
            i, j = j, i
        out[m] = [ids[i], ids[j], current[i, j], sizes[i] + sizes[j]]
        merged = (sizes[i] * current[i] + sizes[j] * current[j]) / (sizes[i] + sizes[j])
        current[i, :] = merged
        current[:, i] = merged
        current[i, i] = np.inf
        current[j, :] = np.inf
        current[:, j] = np.inf
        sizes[i] += sizes[j]
        ids[i] = n + m
        active[j] = False
    return out


def greedy_leaf_order(tree, distances):
    """
    Orders the leaves by flipping every internal node so as to minimize the
    distance between the two rows that become adjacent at its boundary. A
    cheap stand-in for scipy's optimal leaf ordering.
    """
    n = tree.shape[0] + 1
    sequences = {i: [i] for i in range(n)}
    for m in range(n - 1):
        left, right = int(tree[m, 0]), int(tree[m, 1])
        first, second = sequences.pop(left), sequences.pop(right)
        best, best_cost = None, np.inf
        for x in (first, first[::-1]):
            for y in (second, second[::-1]):
                cost = distances[x[-1], y[0]]
                if cost < best_cost:
                    best, best_cost = x + y, cost
        sequences[n + m] = best
    return sequences[2 * n - 2]


def cut_tree(tree, threshold):
    """
    Flat clusters obtained by keeping only the merges below `threshold`.
    Labels are 1-based, as in scipy's fcluster.
    """
    n = tree.shape[0] + 1
    parent = list(range(2 * n - 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for m in range(n - 1):
        if tree[m, 2] <= threshold:
            parent[find(int(tree[m, 0]))] = n + m
            parent[find(int(tree[m, 1]))] = n + m
    seen, labels = {}, np.zeros(n, dtype=int)
    for i in range(n):
        labels[i] = seen.setdefault(find(i), len(seen) + 1)
    return labels


def record_distances(records, args):
    """
    Returns (distances,index,counts): the matrix of the distances between the
    distinct keys of the records, the key of every record, and the number of
    records per key. `--metric profile` has no exact notion of a distinct key,
    so there every record is its own key, and the matrix is quadratic in the
    number of records.
    """
    if args.metric == "edit":
        return edit_distance_matrix(records, args)
    colors, orientations, covered = record_profiles(records, args.bins)
    distances = profile_distance_matrix(
        records,
        colors,
        orientations,
        covered,
        args.w_orientation,
        args.gap_penalty,
        args.len_weight,
    )
    n = len(records)
    return distances.astype(np.float32), np.arange(n), np.ones(n, dtype=int)


def condensed_distances(distances):
    """
    The upper triangle of `distances`, as the double-precision vector scipy's
    linkage expects. Built row by row on purpose: `squareform()` would first
    make a double copy of the whole square matrix, which doubles the peak
    memory of the largest cohorts for nothing.
    """
    n = distances.shape[0]
    out = np.empty(n * (n - 1) // 2, dtype=np.float64)
    position = 0
    for row in range(n - 1):
        size = n - row - 1
        out[position : position + size] = distances[row, row + 1 :]
        position += size
    return out


def cluster_records(distances, args, wants_order=True):
    """
    Returns (order,labels,tree): the row order to draw, the flat cluster label
    of every record, and the tree these come from.
    """
    n = distances.shape[0]
    # The leaf ordering only decides the order of the rows of an image, so
    # there is nothing to gain from it when no image needs it. It is cubic and
    # single-threaded, and dominates everything else past a few thousand keys,
    # so it is also skipped above --olo-max.
    with_ordering = wants_order and n <= args.olo_max
    if HAS_SCIPY:
        condensed = condensed_distances(distances)
        tree = linkage(condensed, method=args.linkage)
        if with_ordering:
            tree = optimal_leaf_ordering(tree, condensed)
        elif wants_order:
            print(
                "%d keys: skipping the optimal leaf ordering (--olo-max is %d)"
                % (n, args.olo_max),
                file=sys.stderr,
            )
        del condensed
        order = list(leaves_list(tree))
        if args.n_clusters > 0:
            labels = fcluster(tree, min(args.n_clusters, n), criterion="maxclust")
        else:
            labels = fcluster(tree, args.cluster_threshold, criterion="distance")
    else:
        tree = average_linkage(distances)
        order = greedy_leaf_order(tree, distances)
        if args.n_clusters > 0:
            k = min(args.n_clusters, n)
            threshold = tree[n - 1 - k, 2] if k < n else -1.0
        else:
            threshold = args.cluster_threshold
        labels = cut_tree(tree, threshold)
    # Relabels clusters top to bottom, so that colors run in drawing order.
    remap, out = {}, np.zeros(n, dtype=int)
    for leaf in order:
        out[leaf] = remap.setdefault(labels[leaf], len(remap) + 1)
    return order, out, tree


# ---------------------------- Dimensionality --------------------------------

def thread_count(threads):
    return threads if threads > 0 else max(1, os.cpu_count() or 1)


def parallel_chunks(total, chunk, work, threads):
    """
    Runs `work(start,stop)` over consecutive chunks of `range(total)`, on
    `threads` threads. Every chunk writes to its own slice of the output, so
    no locking is needed, and numpy releases the interpreter lock inside both
    its elementwise loops and BLAS, so the threads do run in parallel.
    """
    bounds = [(start, min(start + chunk, total)) for start in range(0, total, chunk)]
    if threads <= 1 or len(bounds) == 1:
        for start, stop in bounds:
            work(start, stop)
        return
    with ThreadPoolExecutor(max_workers=min(threads, len(bounds))) as pool:
        list(pool.map(lambda bound: work(*bound), bounds))


def distance_block(first, second):
    """
    Distances between the rows of `first` and the rows of `second`, via the
    Gram matrix rather than by broadcasting: the multiplication goes to BLAS,
    and no temporary of the size of the two clouds is built.
    """
    out = first.dot(second.T)
    out *= -2.0
    out += (first ** 2).sum(1)[:, None]
    out += (second ** 2).sum(1)[None, :]
    np.maximum(out, 0.0, out=out)
    return np.sqrt(out, out=out)


def coordinate_distances(coordinates):
    return distance_block(coordinates, coordinates)


def chunk_size(n, threads):
    """Chunks small enough to fill every thread, large enough to amortize."""
    return int(min(max(n // max(threads * 4, 1), 128), 2048))


def double_centered(distances):
    """
    The matrix whose eigenvectors are the principal coordinates: minus half
    the doubly centered matrix of the squared distances.
    """
    squared = np.asarray(distances, dtype=np.float64) ** 2
    columns = squared.mean(axis=0)
    gram = -0.5 * (
        squared - columns[None, :] - squared.mean(axis=1)[:, None] + squared.mean()
    )
    return (gram + gram.T) / 2, columns


def principal_coordinates(distances, dims=2):
    """
    Classical multidimensional scaling (Torgerson scaling, a.k.a. principal
    coordinates analysis): the projection is read off the eigenvectors of the
    doubly centered matrix of the squared distances. On a precomputed distance
    matrix this is exactly what PCA is: with Euclidean distances it returns
    the principal components of the (unavailable) coordinates, and the
    eigenvalues are the variance each axis carries.

    Returns (coordinates,eigenvalues), the eigenvalues sorted decreasingly.
    """
    gram, _ = double_centered(distances)
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    top = order[:dims]
    return vectors[:, top] * np.sqrt(np.maximum(values[top], 0)), values[order]


def classical_mds(distances, dims=2):
    return principal_coordinates(distances, dims)[0]


def explained_variance(values, dims):
    """
    Fraction of the total variance carried by each of the first `dims` axes.
    Only the positive eigenvalues count: the negative ones are the part of the
    distances that no Euclidean space reproduces.
    """
    total = np.maximum(values, 0).sum()
    if total <= 0:
        return np.zeros(dims)
    return np.maximum(values[:dims], 0) / total


def smacof(distances, start, iterations, threads=1, tolerance=1e-7):
    """
    Metric MDS by majorization (SMACOF): minimizes the raw stress
    sum_ij (|xi-xj| - dij)^2 by repeated Guttman transforms. Unlike classical
    MDS it optimizes the distances themselves, and unlike t-SNE or UMAP it has
    no notion of neighborhood: it is the projection that reproduces the
    distance matrix as closely as 2D allows.

    Every iteration is computed by blocks of rows, in parallel: a block needs
    only its own rows of the Guttman transform, which are then multiplied by
    the whole configuration. The full n x n transform is never materialized.
    """
    n = distances.shape[0]
    # Single precision throughout: every iteration sweeps the whole matrix a
    # few times, so the loop is bound by memory bandwidth rather than by
    # arithmetic, and the extra digits buy nothing in a figure.
    distances = np.asarray(distances, dtype=np.float32)
    coordinates = start.astype(np.float32)
    chunk = chunk_size(n, threads)
    n_chunks = (n + chunk - 1) // chunk
    previous = np.inf
    for _ in range(iterations):
        transformed = np.empty_like(coordinates)
        stresses = np.zeros(n_chunks)

        def work(start_row, stop_row):
            current = distance_block(coordinates[start_row:stop_row], coordinates)
            target = distances[start_row:stop_row]
            stresses[start_row // chunk] = ((current - target) ** 2).sum()
            with np.errstate(divide="ignore", invalid="ignore"):
                block = np.where(current > 1e-12, -target / np.where(current > 0, current, 1), 0.0)
            rows = np.arange(stop_row - start_row)
            block[rows, start_row + rows] = 0.0
            block[rows, start_row + rows] = -block.sum(axis=1)
            transformed[start_row:stop_row] = block.dot(coordinates)

        parallel_chunks(n, chunk, work, threads)
        coordinates = transformed / n
        stress = stresses.sum() / 2
        if previous - stress <= tolerance * max(previous, 1e-12):
            break
        previous = stress
    return coordinates


def landmark_coordinates(distances, n_landmarks, dims, threads=1):
    """
    Principal coordinates of a subset of `n_landmarks` points, with every
    remaining point placed by triangulation from its distances to the
    landmarks alone (de Silva and Tenenbaum). Cost is linear in the number of
    points instead of quadratic, which is what makes tens of thousands of
    points feasible.

    Landmarks are taken by striding the points, which spreads them over the
    data without drawing random numbers.

    Returns (coordinates,eigenvalues,landmark_coordinates,landmarks,anchor).
    """
    n = distances.shape[0]
    landmarks = np.unique(np.linspace(0, n - 1, min(n_landmarks, n)).astype(int))
    anchor = distances[np.ix_(landmarks, landmarks)].astype(np.float64)

    # The triangulation below is consistent with this configuration of the
    # landmarks, and with this one only.
    gram, means = double_centered(anchor)
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    top = order[: min(dims, len(landmarks))]
    scale = np.sqrt(np.maximum(values[top], 1e-12))
    classical = vectors[:, top] * scale
    pseudo = (vectors[:, top] / scale).T  # dims x n_landmarks

    # Every point is triangulated from its distances to the landmarks alone,
    # so the points are independent and the whole loop is parallel.
    out = np.empty((n, len(top)))

    def work(start, stop):
        block = distances[start:stop][:, landmarks].astype(np.float64) ** 2
        block -= means[None, :]
        out[start:stop] = -0.5 * block.dot(pseudo.T)

    parallel_chunks(n, chunk_size(n, threads), work, threads)
    return out, values[order], classical, landmarks, anchor


def landmark_mds(distances, n_landmarks, iterations, threads=1, dims=2):
    """
    Landmark MDS: the landmark configuration of `landmark_coordinates()` is
    refined by majorization, and the triangulated cloud is carried along.
    """
    out, _, classical, landmarks, anchor = landmark_coordinates(
        distances, n_landmarks, dims, threads
    )
    refined = smacof(anchor, classical, iterations, threads).astype(np.float64)
    source, target = classical - classical.mean(0), refined - refined.mean(0)
    left, singular, right = np.linalg.svd(source.T.dot(target))
    rotation = left.dot(right)
    factor = singular.sum() / max((source ** 2).sum(), 1e-12)
    out = (out - classical.mean(0)).dot(rotation) * factor + refined.mean(0)
    out[landmarks] = refined
    return out


def embedding_quality(distances, coordinates, sample):
    """
    Returns (stress1,correlation): Kruskal's stress-1, for which <0.05 is
    excellent, <0.10 good and >0.20 poor; and the correlation between the
    original and the projected distances (the Shepard diagram). Both are
    estimated on `sample` points at most, since they are quadratic.
    """
    n = distances.shape[0]
    if n < 2:
        return 0.0, 1.0
    taken = np.unique(np.linspace(0, n - 1, min(sample, n)).astype(int))
    reduced = distances[np.ix_(taken, taken)].astype(np.float64)
    upper = np.triu_indices_from(reduced, k=1)
    original = reduced[upper]
    projected = coordinate_distances(coordinates[taken])[upper]
    total = (original ** 2).sum()
    stress = np.sqrt(((projected - original) ** 2).sum() / total) if total > 0 else 0.0
    if original.std() == 0 or projected.std() == 0:
        return stress, 1.0
    return stress, float(np.corrcoef(original, projected)[0, 1])


def embed_records(distances, args):
    """
    Projects the distinct keys to 2D. Beyond `--landmarks` points the full
    majorization is replaced by landmark MDS, which is linear rather than
    quadratic in the number of points.
    """
    n = distances.shape[0]
    threads = thread_count(args.threads)
    if args.embed == "pca":
        # Four axes, shown as two panels. No majorization here: majorization
        # optimizes the distances but destroys the ordering and the meaning of
        # the axes, and an axis is exactly what a principal component is.
        if n > args.landmarks:
            coordinates, values = landmark_coordinates(
                distances, args.landmarks, PCA_DIMS, threads
            )[:2]
            method = "pca on %d landmarks, %d threads" % (min(args.landmarks, n), threads)
        else:
            coordinates, values = principal_coordinates(distances, PCA_DIMS)
            method = "pca, %d threads" % threads
        if coordinates.shape[1] < PCA_DIMS:  # Fewer points than axes asked for.
            coordinates = np.pad(
                coordinates, ((0, 0), (0, PCA_DIMS - coordinates.shape[1]))
            )
        fractions = explained_variance(values, PCA_DIMS)
        print(
            "%s of %d points: variance %s"
            % (method, n, " ".join("PC%d %.1f%%" % (i + 1, 100 * f) for i, f in enumerate(fractions))),
            file=sys.stderr,
        )
        stress, correlation = embedding_quality(distances, coordinates[:, :2], args.stress_sample)
        print(
            "  PC1-PC2 alone: stress-1 %.3f, distance correlation %.3f" % (stress, correlation),
            file=sys.stderr,
        )
        return coordinates
    if args.embed == "tsne":
        try:
            from sklearn.manifold import TSNE
        except ImportError:
            sys.exit("--embed tsne needs scikit-learn: pip install scikit-learn")
        perplexity = min(args.perplexity, max(5.0, (n - 1) / 3.0))
        coordinates = TSNE(
            n_components=2,
            metric="precomputed",
            init="random",
            perplexity=perplexity,
            random_state=0,
        ).fit_transform(distances.astype(np.float64))
        method = "tsne"
    elif n > args.landmarks:
        coordinates = landmark_mds(distances, args.landmarks, args.mds_iterations, threads)
        method = "landmark mds on %d landmarks, %d threads" % (min(args.landmarks, n), threads)
    else:
        square = distances.astype(np.float64)
        coordinates = classical_mds(square)
        if args.embed == "mds":  # Refines the classical solution by majorization.
            coordinates = smacof(square, coordinates, args.mds_iterations, threads)
        method = "%s, %d threads" % (args.embed, threads)
    stress, correlation = embedding_quality(distances, coordinates, args.stress_sample)
    print(
        "%s of %d points: stress-1 %.3f, distance correlation %.3f"
        % (method, n, stress, correlation),
        file=sys.stderr,
    )
    return coordinates


# ------------------------------- Drawing ------------------------------------

def uncovered_intervals(segments, alt_length):
    """
    Complement, inside [0,alt_length], of the union of the ALT intervals of
    `segments`. Returns a list of (start,end) in ALT coordinates.
    """
    covered = sorted((segment[4] - 1, segment[5]) for segment in segments)
    merged = []
    for start, end in covered:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    out = []
    cursor = 0
    for start, end in merged:
        if start > cursor:
            out.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < alt_length:
        out.append((cursor, alt_length))
    return out


def draw_record(band, record, colors):
    """
    Rasterizes one record into `band`, an RGB view of the pixels of its row.
    Drawing every record into a single image, instead of adding two matplotlib
    artists per alignment, is what keeps a figure of tens of thousands of rows
    feasible.

    `colors` is the colormap sampled into a lookup table.
    """
    height, width = band.shape[0], band.shape[1]
    alt_length = max(record["alt_length"], 1)
    rows = np.arange(height) - (height - 1) / 2.0
    body_half = max(0.26 * height, 0.5)
    head_half = max(0.44 * height, 0.5)

    # Uncovered portions of ALT: a thin horizontal line.
    line_half = max(height / 24.0, 0.5)
    for start, end in uncovered_intervals(record["segments"], alt_length):
        left = min(int(start / alt_length * width), width - 1)
        right = max(min(int(round(end / alt_length * width)), width), left + 1)
        band[np.abs(rows) <= line_half, left:right] = 140

    for _, seg_start, seg_end, ori, alt_start, alt_end in record["segments"]:
        left = min(int((alt_start - 1) / alt_length * width), width - 1)
        right = max(min(int(round(alt_end / alt_length * width)), width), left + 1)
        columns = np.arange(left, right)
        head = max(min(MAX_HEAD_WIDTH * width, MAX_HEAD_FRACTION * (right - left)), 1.0)
        forward = ori == "+"
        # Distance of every column from the tip of the arrow, in pixels.
        to_tip = (right - columns) if forward else (columns - left + 1)
        half = np.where(to_tip > head, body_half, head_half * to_tip / head)

        # Subgradient of the alignment's reference interval, inside the fixed
        # rainbow gradient of the whole remap_coords interval.
        start_fraction, end_fraction = color_fractions(record, seg_start, seg_end)
        if forward:  # Left-to-right in ALT: same order as in the reference.
            gradient = np.linspace(start_fraction, end_fraction, len(columns))
        else:  # Left-to-right in ALT: reverse order.
            gradient = np.linspace(end_fraction, start_fraction, len(columns))
        rgb = colors[np.clip((gradient * (len(colors) - 1)).astype(int), 0, len(colors) - 1)]

        inside = np.abs(rows)[:, None] <= half[None, :]
        target = band[:, left:right]
        target[inside] = np.broadcast_to(rgb, (height, len(columns), 3))[inside]


def track_figure_height(n_records, args):
    """
    Height of the figure of the tracks, in inches, and height of a row, in
    pixels. Both are capped: at tens of thousands of records the requested
    height would be hundreds of inches, and rendering a canvas that size costs
    minutes and gigabytes.
    """
    height = max(2.0, args.row_height * n_records + 1.0)
    if height > args.max_height:
        height = args.max_height
        print(
            "Figure capped to %g inches: rows are thinner than requested" % height,
            file=sys.stderr,
        )
    row_height = max(int((height - 1.0) * args.dpi) // max(n_records, 1), 1)
    return height, row_height


def render_tracks(records, order, row_height, args, cmap):
    """
    Rasterizes every record of `order` into one RGB image, one row per record,
    every row of the same width. A single image, instead of matplotlib artists
    per alignment, is what makes tens of thousands of rows feasible.
    """
    width = max(int(round(args.width * args.dpi)), 16)
    budget = int(args.max_pixels * 1e6)
    if len(order) * row_height * width > budget:  # Keeps the image within memory.
        row_height = max(budget // (max(len(order), 1) * width), 1)
        print(
            "Rows reduced to %d pixel(s) to stay under --max-pixels" % row_height,
            file=sys.stderr,
        )
    out = np.full((len(order) * row_height, width, 3), 255, dtype=np.uint8)
    colors = np.clip(cmap(np.linspace(0, 1, 512))[:, :3] * 255, 0, 255).astype(np.uint8)
    for position, leaf in enumerate(order):
        top = position * row_height
        draw_record(out[top : top + row_height], records[leaf], colors)
    return out


def source_profile(record, n_bins):
    """
    Coverage of the `remap_coords` interval, normalized to [0,1] and quantized
    uniformly into `n_bins`: every alignment adds one unit of mass to every
    quantum its reference interval overlaps. Where the tracks show what
    happens along ALT, this shows which parts of the source interval are used,
    and how many times.
    """
    _, ref_start, ref_end = record["coords"]
    span = max(ref_end - ref_start, 1)
    out = np.zeros(n_bins, dtype=np.float32)
    for _, seg_start, seg_end, _, _, _ in record["segments"]:
        low = (seg_start - ref_start) / span
        high = (seg_end - ref_start) / span
        if high < low:
            low, high = high, low
        first = int(np.clip(np.floor(low * n_bins), 0, n_bins - 1))
        last = int(np.clip(np.ceil(high * n_bins) - 1, 0, n_bins - 1))
        out[first : last + 1] += 1.0
    return out


def source_profiles(records, n_bins):
    """
    One probability distribution over the quantized source interval per
    record: the coverage of `source_profile()`, divided by its total mass.
    """
    out = np.empty((len(records), n_bins), dtype=np.float32)
    for index, record in enumerate(records):
        out[index] = source_profile(record, n_bins)
    totals = out.sum(axis=1, keepdims=True)
    np.divide(out, np.maximum(totals, 1e-12), out=out)
    return out


def euclidean_distances(points, threads):
    """
    Matrix of the Euclidean distances between the rows of `points`, computed
    by blocks of rows in parallel.
    """
    n = points.shape[0]
    out = np.empty((n, n), dtype=np.float32)

    def work(start, stop):
        out[start:stop] = distance_block(points[start:stop], points)

    parallel_chunks(n, KEY_CHUNK, work, threads)
    return out


def render_source(profiles, order, row_height, width, cmap, vmax):
    """
    Rasterizes the distributions of `order`, one row each, coloring a quantum
    by its probability in multiples of the uniform distribution, so that rows
    of different peakedness stay comparable: 1 is a flat distribution, and
    everything at or above `vmax` saturates.
    """
    n_bins = profiles.shape[1]
    repeat = max(int(round(width / n_bins)), 1)
    colors = np.clip(cmap(np.linspace(0, 1, 512))[:, :3] * 255, 0, 255).astype(np.uint8)
    scaled = np.clip(profiles[order] * n_bins / max(vmax, 1e-9), 0.0, 1.0)
    indices = (scaled * (len(colors) - 1)).astype(np.int32)
    image = colors[indices]  # rows x n_bins x 3
    return np.repeat(np.repeat(image, row_height, axis=0), repeat, axis=1)


def draw_source(path, records, args, palette):
    """
    One row per call, showing how the alignments cover the source interval,
    with the rows ordered by a hierarchical clustering of the distributions
    under the Euclidean distance.
    """
    threads = thread_count(args.threads)
    profiles = source_profiles(records, args.source_bins)

    # Identical distributions are one point for the clustering, exactly as the
    # identical strings are for the edit distance.
    unique, key_of = {}, []
    for row in profiles:
        key_of.append(unique.setdefault(row.tobytes(), len(unique)))
    key_of = np.array(key_of)
    members = [[] for _ in range(len(unique))]
    for record, key in enumerate(key_of):
        members[key].append(record)
    representatives = np.array([group[0] for group in members])
    print(
        "%d records, %d distinct distributions over %d quanta"
        % (len(records), len(unique), args.source_bins),
        file=sys.stderr,
    )

    if len(unique) >= 3 and not args.no_cluster:
        distances = euclidean_distances(profiles[representatives], threads)
        key_order, key_labels, tree = cluster_records(distances, args)
        # The tree is on the Euclidean distance between distributions, whose
        # scale has nothing in common with the edit distance between strings,
        # so --cluster-threshold does not apply: the cut is taken relative to
        # the tallest merge instead.
        threshold = args.source_threshold
        if threshold <= 0:
            threshold = SOURCE_THRESHOLD_FRACTION * float(tree[:, 2].max())
        flat = fcluster(tree, threshold, criterion="distance") if HAS_SCIPY else cut_tree(
            tree, threshold
        )
        remap, key_labels = {}, np.zeros(len(unique), dtype=int)
        for leaf in key_order:  # Colors run in drawing order.
            key_labels[leaf] = remap.setdefault(flat[leaf], len(remap) + 1)
        sizes = np.bincount(key_labels)[1:]
        print(
            "%d clusters of distributions, sizes %s"
            % (len(sizes), " ".join(str(size) for size in sizes[:200])),
            file=sys.stderr,
        )
    else:
        key_order = list(range(len(unique)))
        key_labels, tree = np.ones(len(unique), dtype=int), None
    order = [record for key in key_order for record in members[key]]

    figure_height, row_height = track_figure_height(len(records), args)
    width = max(int(round(args.width * args.dpi)), 16)
    budget = int(args.max_pixels * 1e6)
    if len(order) * row_height * width > budget:
        row_height = max(budget // (max(len(order), 1) * width), 1)
    raster = render_source(
        profiles, order, row_height, width, plt.get_cmap(args.source_cmap), args.source_vmax
    )

    if tree is None or not args.dendrogram:
        plt.imsave(path, raster)
    else:
        figure = plt.figure(
            figsize=(args.width, figure_height), dpi=args.dpi, layout="constrained"
        )
        grid = figure.add_gridspec(1, 2, width_ratios=[0.13, 1], wspace=0.01)
        ax_tree, ax = figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])
        ax.imshow(raster, extent=(0.0, 1.0, float(len(order)), 0.0), aspect="auto")
        ax.set_xlabel("Position in the source interval (normalized)")
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.set_yticks([])
        ax.set_ylim(len(order), 0.0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        # A leaf sits at the middle of the block of rows of its distribution.
        key_row, position = {}, 0
        for key in key_order:
            key_row[key] = position + len(members[key]) / 2.0
            position += len(members[key])
        maximum = draw_dendrogram(ax_tree, tree, key_row, key_labels, palette)
        ax_tree.set_xlim(maximum * 1.05, 0)
        ax_tree.set_ylim(len(order), 0.0)
        ax_tree.set_yticks([])
        ax_tree.tick_params(axis="x", labelsize=6)
        ax_tree.set_xlabel("Distance", fontsize=7)
        for side in ("top", "right", "left"):
            ax_tree.spines[side].set_visible(False)
        figure.savefig(path, dpi=args.dpi, bbox_inches="tight")
        plt.close(figure)
    print(
        "Saved %d source distributions to %s (color saturates at %g times uniform)"
        % (len(records), path, args.source_vmax),
        file=sys.stderr,
    )


def draw_dendrogram(ax, tree, row_of, labels, palette):
    """
    Draws `tree` with its leaves aligned to the rows of the main axes: leaf `i`
    is drawn at height `row_of[i]`. Branches inside a flat cluster take the
    color of the cluster.
    """
    n = tree.shape[0] + 1
    coordinates = {i: (0.0, float(row_of[i])) for i in range(n)}
    cluster_of = {i: int(labels[i]) for i in range(n)}
    for m in range(n - 1):
        left, right = int(tree[m, 0]), int(tree[m, 1])
        (x_left, y_left), (x_right, y_right) = coordinates[left], coordinates[right]
        distance = float(tree[m, 2])
        same = cluster_of.get(left) is not None and cluster_of.get(left) == cluster_of.get(right)
        color = palette[(cluster_of[left] - 1) % len(palette)] if same else "0.45"
        ax.plot(
            [x_left, distance, distance, x_right],
            [y_left, y_left, y_right, y_right],
            color=color,
            linewidth=0.6,
            solid_joinstyle="miter",
        )
        coordinates[n + m] = (distance, (y_left + y_right) / 2)
        cluster_of[n + m] = cluster_of[left] if same else None
    return float(tree[:, 2].max())


def union_length(intervals):
    """
    Total length of the union of a set of intervals, which is what makes the
    length-weighted statistics well defined when alignments overlap in ALT.
    """
    total, current_start, current_end = 0, None, None
    for start, end in sorted(intervals):
        if current_end is None or start > current_end:
            if current_end is not None:
                total += current_end - current_start
            current_start, current_end = start, end
        else:
            current_end = max(current_end, end)
    if current_end is not None:
        total += current_end - current_start
    return total


def covered_fraction(record):
    """
    Fraction of ALT that the union of the alignments covers, in [0,1].
    """
    alt_length = max(record["alt_length"], 1)
    uncovered = sum(
        end - start for start, end in uncovered_intervals(record["segments"], alt_length)
    )
    return 1.0 - uncovered / alt_length


def orientation_fractions(record):
    """
    Fraction of the whole ALT covered by the forward alignments, and by the
    reverse-complemented ones. Both are unions, so overlapping alignments of
    the same orientation are not counted twice; the two can still overlap each
    other, so the pair does not have to sum to the covered fraction.
    """
    alt_length = max(record["alt_length"], 1)
    by_strand = []
    for wanted in ("+", "-"):
        by_strand.append(
            union_length(
                [
                    (segment[4] - 1, segment[5])
                    for segment in record["segments"]
                    if segment[3] == wanted
                ]
            )
            / alt_length
        )
    return min(by_strand[0], 1.0), min(by_strand[1], 1.0)


def reverse_fraction(record):
    """
    Fraction of the covered part of ALT that comes from reverse-complemented
    alignments, weighted by length rather than counted per alignment: one long
    inversion and twenty short ones are very different events, and counting
    alignments cannot tell them apart.
    """
    intervals = [(segment[4] - 1, segment[5]) for segment in record["segments"]]
    covered = union_length(intervals)
    if covered <= 0:
        return 0.0
    reverse = union_length(
        [
            (segment[4] - 1, segment[5])
            for segment in record["segments"]
            if segment[3] == "-"
        ]
    )
    return min(reverse / covered, 1.0)


def orientation_counts(records):
    """
    Number of alignments in each orientation, per record.
    """
    forward = np.array(
        [sum(1 for segment in record["segments"] if segment[3] == "+") for record in records]
    )
    total = np.array([len(record["segments"]) for record in records])
    return forward, total - forward


def bare_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def capped_ticks(ax, cap):
    """
    Integer ticks up to `cap`, the last one marked as an overflow.
    """
    ticks = np.arange(0, cap + 1, max(1, cap // 10))
    labels = [str(tick) for tick in ticks]
    if ticks[-1] == cap:  # The last row and column absorb everything above.
        labels[-1] = "%d+" % cap
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)


def draw_orientation_heatmap(figure, ax, forward, reverse, cap):
    """
    Joint distribution of the number of alignments per orientation.
    """
    grid = np.zeros((cap + 1, cap + 1), dtype=np.int64)
    np.add.at(grid, (np.minimum(reverse, cap), np.minimum(forward, cap)), 1)
    image = ax.imshow(
        np.ma.masked_equal(grid, 0),  # Empty cells stay white.
        origin="lower",
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(grid.max(), 2)),
        aspect="equal",
        interpolation="nearest",
        extent=(-0.5, cap + 0.5, -0.5, cap + 0.5),
    )
    ax.set_xlabel("Number of FWD alignments")
    ax.set_ylabel("Number of RC alignments")
    capped_ticks(ax, cap)
    figure.colorbar(image, ax=ax, label="Calls", shrink=0.85)


def draw_coverage_heatmap(figure, ax, forward, reverse, bins):
    """
    Joint distribution of the fraction of ALT covered by each orientation.
    The same picture as the heatmap of the counts, but weighted by how much
    sequence each alignment actually explains.
    """
    edges = np.linspace(0.0, 1.0, bins + 1)
    grid, _, _ = np.histogram2d(forward, reverse, bins=[edges, edges])
    image = ax.pcolormesh(
        edges,
        edges,
        np.ma.masked_equal(grid.T, 0),
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(grid.max(), 2)),
    )
    # Below this line ALT is not fully covered; on it, the two orientations
    # tile ALT exactly; above it, they overlap each other.
    ax.plot([0, 1], [1, 0], color="0.4", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Fraction of ALT covered by FWD alignments")
    ax.set_ylabel("Fraction of ALT covered by RC alignments")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    figure.colorbar(image, ax=ax, label="Calls", shrink=0.85)


def draw_length_heatmap(figure, ax, lengths, counts, cap, bins, max_length):
    """
    Number of alignments against ALT length: tells whether a call is shredded
    because it is long, or because of where it is. Calls longer than
    `max_length` are left out of this panel rather than piled into its last
    bin, which would put a spike where the distribution has a tail; how many
    were dropped is reported.
    """
    if max_length > 0:
        inside = lengths <= max_length
        dropped = int((~inside).sum())
        if dropped:
            print(
                "  %d call(s) longer than %d bp are not in the length panel"
                % (dropped, max_length),
                file=sys.stderr,
            )
        lengths, counts = lengths[inside], counts[inside]
        if lengths.size == 0:
            ax.set_axis_off()
            return
    low, high = float(lengths.min()), float(lengths.max())
    if high - low < 1e-9:  # Every call of the same length.
        low, high = low - 0.5, high + 0.5
    x_edges = np.linspace(low, high, bins + 1)
    y_edges = np.arange(-0.5, cap + 1.5)
    grid, _, _ = np.histogram2d(lengths, np.minimum(counts, cap), bins=[x_edges, y_edges])
    image = ax.pcolormesh(
        x_edges,
        y_edges,
        np.ma.masked_equal(grid.T, 0),
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(grid.max(), 2)),
    )
    ax.set_xlabel("ALT length (bp)")
    ax.set_ylabel("Number of alignments")
    ax.set_xlim(low, high)
    ax.set_ylim(-0.5, cap + 0.5)
    figure.colorbar(image, ax=ax, label="Calls", shrink=0.85)


def alignment_alt_lengths(records):
    """
    Length in ALT of every alignment of every record, one entry per alignment
    rather than one per call: the other panels say how ALT is partitioned, this
    one says how long the pieces of the partition are. Returns the lengths of
    the forward alignments and those of the reverse-complemented ones,
    separately.
    """
    out = {"+": [], "-": []}
    for record in records:
        for segment in record["segments"]:
            out.setdefault(segment[3], []).append(segment[5] - segment[4] + 1)
    return (
        np.array(out["+"], dtype=np.float64),
        np.array(out["-"], dtype=np.float64),
    )


def step_histograms(ax, forward_lengths, reverse_lengths, edges, linewidth):
    """
    The two orientations as two outlines over the same bins, so that their
    heights can be compared bin by bin. Outlines rather than filled bars: two
    translucent fills would overlap into a third color that reads as a third
    series.
    """
    for lengths, color, label in (
        (forward_lengths, FORWARD_COLOR, "FWD"),
        (reverse_lengths, REVERSE_COLOR, "RC"),
    ):
        if lengths.size == 0:  # No entry in the legend for an empty series.
            continue
        ax.hist(
            lengths,
            bins=edges,
            histtype="step",
            linewidth=linewidth,
            color=color,
            label=label,
        )


def draw_alignment_length_histogram(ax, forward_lengths, reverse_lengths, bins):
    """
    Distribution of the length in ALT of a single alignment, one histogram per
    orientation. The lengths are binned linearly; only the counts are on a log
    axis, since the bulk would otherwise hide the whole tail.

    An inset repeats the two histograms over [0,INSET_MAX_LENGTH], with the
    same number of bins over that shorter range: on the full range a few very
    long alignments stretch the axis and squeeze the bulk of the distribution
    into its first bins.
    """
    both = np.concatenate([forward_lengths, reverse_lengths])
    if both.size == 0:
        ax.set_axis_off()
        return
    low, high = float(both.min()), float(both.max())
    if high - low < 1e-9:  # Every alignment of the same length.
        low, high = low - 0.5, high + 0.5
    step_histograms(ax, forward_lengths, reverse_lengths, np.linspace(low, high, bins + 1), 1.4)
    ax.set_yscale("log")
    ax.set_xlabel("ALT length of an alignment (bp)")
    ax.set_ylabel("Alignments (log scale)")
    ax.set_xlim(low, high)
    # The legend moves out of the way of the inset, into the corner that the
    # decay of the distribution leaves empty.
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    bare_axes(ax)

    # Flush with the top-right corner of the panel: width and height first,
    # then the origin that puts the far corner at (1,1).
    width, height = 0.58, 0.55
    inset = ax.inset_axes([1.0 - width, 1.0 - height, width, height])
    step_histograms(
        inset,
        forward_lengths,
        reverse_lengths,
        np.linspace(0.0, float(INSET_MAX_LENGTH), bins + 1),
        1.0,
    )
    inset.set_yscale("log")
    inset.set_xlim(0.0, float(INSET_MAX_LENGTH))
    # The caption goes inside: a title would sit above the top spine, outside
    # the panel that the inset is flush with.
    inset.text(
        0.02,
        0.97,
        "0-%d bp, %d bins" % (INSET_MAX_LENGTH, bins),
        transform=inset.transAxes,
        ha="left",
        va="top",
        fontsize=7,
    )
    inset.tick_params(labelsize=6)
    # An opaque, fully framed box: without it the curves of the inset and
    # those of the panel underneath read as one tangle.
    inset.set_facecolor("white")
    for side in ("top", "right", "bottom", "left"):
        inset.spines[side].set_visible(True)
        inset.spines[side].set_color("0.6")
        inset.spines[side].set_linewidth(0.8)


def draw_summary(path, records, args):
    """
    Six marginal summaries of the whole cohort, in one figure: how much of
    ALT the alignments cover, the number of alignments per orientation, the
    fraction of ALT covered per orientation, the number of alignments against
    ALT length, how much of ALT is inverted by length rather than by count,
    and the ALT length of a single alignment, per orientation.
    """
    fractions = np.array([covered_fraction(record) for record in records])
    inverted = np.array([reverse_fraction(record) for record in records])
    lengths = np.array([record["alt_length"] for record in records])
    covered = np.array([orientation_fractions(record) for record in records])
    forward_lengths, reverse_lengths = alignment_alt_lengths(records)
    forward, reverse = orientation_counts(records)
    counts = forward + reverse

    cap = args.max_alignments
    if cap <= 0:  # Keeps the heatmaps readable when a few calls are shredded.
        cap = int(np.clip(np.percentile(counts, 99.5), 1, 40))
    orientation_cap = max(
        int(np.clip(np.percentile(np.concatenate([forward, reverse]), 99.5), 1, 40)), 1
    )

    figure, axes = plt.subplots(
        3,
        2,
        figsize=(2.1 * args.scatter_size, 2.6 * args.scatter_size),
        dpi=args.dpi,
        layout="constrained",
    )
    axes[0, 0].hist(fractions, bins=args.summary_bins, range=(0.0, 1.0), color="0.35")
    axes[0, 0].set_yscale("log")  # The bulk would otherwise hide the whole tail.
    axes[0, 0].set_xlabel("Fraction of ALT covered by alignments")
    axes[0, 0].set_ylabel("Calls (log scale)")
    axes[0, 0].set_xlim(0.0, 1.0)
    bare_axes(axes[0, 0])

    draw_orientation_heatmap(figure, axes[0, 1], forward, reverse, orientation_cap)
    draw_coverage_heatmap(figure, axes[1, 0], covered[:, 0], covered[:, 1], args.summary_bins)
    draw_length_heatmap(
        figure, axes[1, 1], lengths, counts, cap, args.summary_bins, args.max_alt_length
    )
    # The panel is drawn on a linear axis, but the correlation stays on the
    # log of the length: lengths span decades, and on a linear scale a handful
    # of very long calls would decide the coefficient on their own.
    logs = np.log10(np.maximum(lengths, 1))

    axes[2, 0].hist(inverted, bins=args.summary_bins, range=(0.0, 1.0), color="0.35")
    axes[2, 0].set_yscale("log")
    axes[2, 0].set_xlabel("Fraction of the covered ALT that is reverse-complemented")
    axes[2, 0].set_ylabel("Calls (log scale)")
    axes[2, 0].set_xlim(0.0, 1.0)
    bare_axes(axes[2, 0])

    draw_alignment_length_histogram(
        axes[2, 1], forward_lengths, reverse_lengths, args.summary_bins
    )

    figure.savefig(path, dpi=args.dpi, bbox_inches="tight")
    plt.close(figure)
    report_summary(
        records,
        fractions,
        inverted,
        lengths,
        logs,
        forward,
        reverse,
        counts,
        covered,
        forward_lengths,
        reverse_lengths,
        path,
    )


def report_summary(
    records,
    fractions,
    inverted,
    lengths,
    logs,
    forward,
    reverse,
    counts,
    covered,
    forward_lengths,
    reverse_lengths,
    path,
):
    """
    The headline numbers behind the four panels.
    """
    total = len(records)
    print(
        "Summary of %d calls: ALT covered median %.3f, mean %.3f; fully covered %.1f%%, "
        "less than half covered %.1f%%"
        % (
            total,
            np.median(fractions),
            fractions.mean(),
            100 * (fractions > 0.99).mean(),
            100 * (fractions < 0.5).mean(),
        ),
        file=sys.stderr,
    )
    print(
        "  alignments per call: median %d, mean %.1f, max %d; one alignment only %.1f%%"
        % (np.median(counts), counts.mean(), counts.max(), 100 * (counts == 1).mean()),
        file=sys.stderr,
    )
    print(
        "  no reverse alignment %.1f%%, no forward alignment %.1f%%, both orientations %.1f%%"
        % (
            100 * (reverse == 0).mean(),
            100 * (forward == 0).mean(),
            100 * ((forward > 0) & (reverse > 0)).mean(),
        ),
        file=sys.stderr,
    )
    if total > 2 and logs.std() > 0 and counts.std() > 0:
        # Does shredding scale with size, or is it a property of the locus?
        correlation = float(np.corrcoef(logs, counts)[0, 1])
        low, high = np.percentile(lengths, 10), np.percentile(lengths, 90)
        short, long_ = counts[lengths <= low], counts[lengths >= high]
        print(
            "  ALT length vs alignments: correlation %.3f on the log of the length; median "
            "%d alignments in the shortest decile (<=%d bp), %d in the longest (>=%d bp)"
            % (correlation, np.median(short), low, np.median(long_), high),
            file=sys.stderr,
        )
    print(
        "  of the whole ALT: forward covers median %.3f, reverse covers median %.3f; the two "
        "orientations overlap in %.1f%% of the calls"
        % (
            np.median(covered[:, 0]),
            np.median(covered[:, 1]),
            100 * (covered.sum(axis=1) > fractions + 0.001).mean(),
        ),
        file=sys.stderr,
    )
    print(
        "  reverse-complemented by length: median %.3f, mean %.3f; none %.1f%%, "
        "more than half %.1f%%, all %.1f%%"
        % (
            np.median(inverted),
            inverted.mean(),
            100 * (inverted <= 0.001).mean(),
            100 * (inverted > 0.5).mean(),
            100 * (inverted >= 0.999).mean(),
        ),
        file=sys.stderr,
    )
    for lengths_of_strand, name in ((forward_lengths, "FWD"), (reverse_lengths, "RC")):
        if lengths_of_strand.size == 0:
            continue
        print(
            "  %d %s alignments, ALT length median %d bp, mean %.1f bp, min %d, max %d; "
            "shorter than 50 bp %.1f%%, longer than 1 kbp %.1f%%"
            % (
                lengths_of_strand.size,
                name,
                np.median(lengths_of_strand),
                lengths_of_strand.mean(),
                lengths_of_strand.min(),
                lengths_of_strand.max(),
                100 * (lengths_of_strand < 50).mean(),
                100 * (lengths_of_strand > 1000).mean(),
            ),
            file=sys.stderr,
        )
    print("Saved the summary of %d calls to %s" % (total, path), file=sys.stderr)


def report_timings():
    total = sum(seconds for _, seconds in TIMINGS)
    print(
        "Timing: %s | total %.1fs"
        % (" ".join("%s %.1fs" % (name, seconds) for name, seconds in TIMINGS), total),
        file=sys.stderr,
    )


def draw_scatter(path, coordinates, counts, labels, palette, clustered, args, pairs=((0, 1),)):
    """
    Writes the projection: one point per distinct key, colored by flat
    cluster, with an area proportional to the number of records it stands for.
    One panel per pair of axes of `pairs`, side by side in one figure.
    """
    figure, axes = plt.subplots(
        1,
        len(pairs),
        figsize=(args.scatter_size * len(pairs), args.scatter_size),
        dpi=args.dpi,
        layout="constrained",
    )
    colors = (
        [palette[(label - 1) % len(palette)] for label in labels]
        if clustered
        else ["0.35"] * len(counts)
    )
    for ax, (first, second) in zip(np.atleast_1d(axes), pairs):
        ax.scatter(
            coordinates[:, first],
            coordinates[:, second],
            s=args.point_size * np.sqrt(counts),
            c=colors,
            alpha=0.85,
            linewidths=0.3,
            edgecolors="0.2",
        )
        ax.set_aspect("equal")  # The axes must not be scaled independently.
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)
        if len(pairs) > 1:  # Only to tell the panels apart.
            ax.set_title("PC%d - PC%d" % (first + 1, second + 1), fontsize=9)
    figure.savefig(path, dpi=args.dpi, bbox_inches="tight")
    plt.close(figure)
    print("Saved %d points to %s" % (len(counts), path), file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Plots truvari anno remap alignments of every VCF in a directory."
    )
    parser.add_argument("input_dir", help="Directory containing .vcf / .vcf.gz files.")
    parser.add_argument("-o", "--output", default="remap_segments.png", help="Output image.")
    parser.add_argument("--cmap", default="gist_rainbow", help="Rainbow colormap to use.")
    parser.add_argument("--row-height", type=float, default=0.30, help="Inches per record.")
    parser.add_argument("--width", type=float, default=14.0, help="Figure width, in inches.")
    parser.add_argument("--dpi", type=int, default=200, help="Output resolution.")
    parser.add_argument(
        "--max-records", type=int, default=0, help="Plots at most this many records (0=all)."
    )
    parser.add_argument(
        "--no-cluster",
        action="store_true",
        help="Keeps the records in file order, instead of clustering them.",
    )
    parser.add_argument(
        "--dendrogram", action="store_true", help="Draws the dendrogram of the clustering."
    )
    parser.add_argument(
        "--labels", action="store_true", help="Writes the origin of every record next to its row."
    )
    parser.add_argument(
        "--cluster-bar",
        action="store_true",
        help="Draws a color bar of the flat clusters, left of the tracks.",
    )
    parser.add_argument(
        "--scatter",
        default=None,
        metavar="FILE",
        help="Also writes a 2D projection of the records to FILE, one point per record, "
        "colored by cluster.",
    )
    parser.add_argument(
        "--embed",
        default="mds",
        choices=("mds", "classical", "pca", "tsne"),
        help="Projection for --scatter: metric MDS by majorization, which minimizes the "
        "error on the distances themselves (default); classical MDS alone; pca, which is "
        "classical MDS on 4 axes, drawn as two panels, PC1-PC2 and PC3-PC4; or t-SNE, "
        "which preserves neighborhoods rather than distances (needs scikit-learn).",
    )
    parser.add_argument(
        "--mds-iterations", type=int, default=300, help="Majorization steps of --embed mds."
    )
    parser.add_argument(
        "--perplexity", type=float, default=30.0, help="Perplexity of --embed tsne."
    )
    parser.add_argument(
        "--scatter-size", type=float, default=7.0, help="Side of the scatter plot, in inches."
    )
    parser.add_argument("--point-size", type=float, default=22.0, help="Area of a single point.")
    parser.add_argument(
        "--no-tracks", action="store_true", help="Skips the figure of the tracks."
    )
    parser.add_argument(
        "--summary",
        default=None,
        metavar="FILE",
        help="Also writes to FILE a figure with six marginal summaries of the cohort: the "
        "fraction of ALT covered by alignments, a heatmap of the number of forward against "
        "reverse-complement alignments, a heatmap of the fraction of ALT that each "
        "orientation covers, a heatmap of the number of alignments against ALT length, "
        "the fraction of ALT that is reverse-complemented, by length, and the distribution "
        "of the ALT length of a single alignment, one histogram per orientation.",
    )
    parser.add_argument(
        "--source",
        default=None,
        metavar="FILE",
        help="Also writes to FILE one row per call, showing how its alignments cover the "
        "remap_coords interval: the interval is normalized to [0,1] and quantized, every "
        "alignment adds one unit of mass to every quantum it overlaps, and the result is "
        "normalized into a probability distribution. Rows are ordered by a hierarchical "
        "clustering of the distributions under the Euclidean distance.",
    )
    parser.add_argument(
        "--source-bins",
        type=int,
        default=256,
        help="Quanta the source interval is partitioned into, for --source.",
    )
    parser.add_argument(
        "--source-cmap", default="magma", help="Colormap of --source."
    )
    parser.add_argument(
        "--source-threshold",
        type=float,
        default=0.0,
        help="Distance at which the tree of --source is cut into flat clusters (0 = %g "
        "times its tallest merge). --cluster-threshold does not apply there: it is on the "
        "edit distance between strings, not on the Euclidean distance between "
        "distributions." % SOURCE_THRESHOLD_FRACTION,
    )
    parser.add_argument(
        "--source-vmax",
        type=float,
        default=4.0,
        help="Probability at which the color of --source saturates, in multiples of the "
        "uniform distribution (1 is a flat distribution over the whole interval).",
    )
    parser.add_argument(
        "--summary-bins",
        type=int,
        default=50,
        help="Bins of the histograms and of the length axis of --summary.",
    )
    parser.add_argument(
        "--max-alignments",
        type=int,
        default=0,
        help="Largest number of alignments shown by the heatmaps of --summary; calls above "
        "it fall in the last row or column (0 = the 99.5th percentile, capped at 40).",
    )
    parser.add_argument(
        "--max-alt-length",
        type=int,
        default=100000,
        help="Largest ALT length shown by the length panel of --summary; longer calls are "
        "left out of that panel, and counted on stderr (0 = no limit).",
    )
    parser.add_argument(
        "--metric",
        default="edit",
        choices=("edit", "profile"),
        help="Row similarity: edit distance between the strings of quantized alignments and "
        "gaps (default), or per-bin disagreement between the rasterized rows.",
    )
    parser.add_argument(
        "--ref-splits",
        type=int,
        default=20,
        help="Splitpoints partitioning remap_coords, for --metric edit. Must stay coarser "
        "than the variation of the breakpoints between records that should cluster "
        "together, unless --sub-cost is graded.",
    )
    parser.add_argument(
        "--alt-chunks",
        type=int,
        default=20,
        help="Chunks partitioning ALT, quantizing the length of a gap, for --metric edit.",
    )
    parser.add_argument(
        "--sub-cost",
        default="unit",
        choices=("unit", "graded"),
        help="Replacement cost in the edit distance: 1 for every distinct pair of characters "
        "(default), or proportional to how far apart the two characters are.",
    )
    parser.add_argument(
        "--bins", type=int, default=128, help="Bins along ALT, for --metric profile."
    )
    parser.add_argument(
        "--w-orientation",
        type=float,
        default=0.4,
        help="Weight of a strand disagreement in the distance; 1-w is the weight of the "
        "difference between reference positions. Used by --metric profile, and by "
        "--metric edit when --sub-cost is graded.",
    )
    parser.add_argument(
        "--gap-penalty",
        type=float,
        default=0.6,
        help="Distance charged to a bin covered in one record and uncovered in the other, "
        "for --metric profile.",
    )
    parser.add_argument(
        "--len-weight",
        type=float,
        default=0.1,
        help="Weight of the ALT-length ratio in the distance (0 to ignore length entirely).",
    )
    parser.add_argument("--linkage", default="average", help="Linkage method (scipy only).")
    parser.add_argument(
        "--olo-max",
        type=int,
        default=1000,
        help="Above this many distinct keys, the optimal leaf ordering, which is cubic and "
        "single-threaded, is replaced by the plain order of the leaves of the tree. It "
        "costs ~1.5s at 1000 keys, ~20s at 2000 and ~90s at 3000.",
    )
    parser.add_argument(
        "--landmarks",
        type=int,
        default=1500,
        help="Above this many distinct keys, --scatter switches to landmark MDS, which is "
        "linear rather than quadratic in the number of points.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=0,
        help="Threads used by the distance matrix and by the projection (0 = every core).",
    )
    parser.add_argument(
        "--stress-sample",
        type=int,
        default=2000,
        help="Points on which the stress of the projection is estimated.",
    )
    parser.add_argument(
        "--max-pixels",
        type=float,
        default=80.0,
        help="Millions of pixels of the raster of the tracks. Rows are made thinner when "
        "the figure would exceed it.",
    )
    parser.add_argument(
        "--max-height",
        type=float,
        default=100.0,
        help="Largest height of the figure of the tracks, in inches. Rows are made thinner "
        "when --row-height times the number of records would exceed it.",
    )
    parser.add_argument(
        "--cluster-threshold",
        type=float,
        default=None,
        help="Distance at which the tree is cut into flat clusters. Defaults to %g for "
        "--metric edit and to %g for --metric profile." % (EDIT_THRESHOLD, PROFILE_THRESHOLD),
    )
    parser.add_argument(
        "--n-clusters", type=int, default=0, help="Cuts into this many clusters instead (0=off)."
    )
    args = parser.parse_args()
    if args.cluster_threshold is None:
        args.cluster_threshold = EDIT_THRESHOLD if args.metric == "edit" else PROFILE_THRESHOLD

    files = sorted(
        set(glob.glob(os.path.join(args.input_dir, "*.vcf")))
        | set(glob.glob(os.path.join(args.input_dir, "*.vcf.gz")))
        | set(glob.glob(os.path.join(args.input_dir, "*.vcf.bgz")))
    )
    if not files:
        sys.exit("No VCF found in %s" % args.input_dir)

    records = []
    with phase("parse"):
        for path in files:
            found = load_records(path)
            print(
                "%s: %d annotated records" % (os.path.basename(path), len(found)),
                file=sys.stderr,
            )
            records.extend(found)
    if not records:
        sys.exit("No record with remap_coords and remap_segments found.")
    if args.max_records > 0 and len(records) > args.max_records:
        print(
            "Plotting the first %d records out of %d" % (args.max_records, len(records)),
            file=sys.stderr,
        )
        records = records[: args.max_records]

    n_records = len(records)
    palette = [plt.get_cmap("tab20")(i / 20) for i in range(20)]
    if args.summary:
        with phase("summary"):
            draw_summary(args.summary, records, args)
    if args.source:
        with phase("source"):
            draw_source(args.source, records, args, palette)

    # The clustering serves the order of the rows and the colors of the
    # scatter plot; with neither figure asked for, the whole distance matrix
    # would be computed for nothing.
    needed = (not args.no_tracks) or args.scatter
    clustered = (not args.no_cluster) and n_records >= 3 and needed
    if not needed:
        report_timings()
        return
    if clustered or args.scatter:
        with phase("distances"):
            distances, key_of, counts = record_distances(records, args)
    else:
        distances, key_of, counts = None, np.arange(n_records), np.ones(n_records, dtype=int)

    # Everything is clustered and projected at the level of the distinct keys;
    # the records of a key are then expanded back, contiguously.
    members = [[] for _ in range(len(counts))]
    for record, key in enumerate(key_of):
        members[key].append(record)
    if clustered:
        if not HAS_SCIPY:
            print("scipy not found: using the built-in clustering.", file=sys.stderr)
        with phase("cluster"):
            key_order, key_labels, tree = cluster_records(distances, args, not args.no_tracks)
        order = [record for key in key_order for record in members[key]]
        labels = key_labels[key_of]
        sizes = np.bincount(labels)[1:]
        print(
            "%d clusters, sizes %s"
            % (len(sizes), " ".join(str(size) for size in sizes[: args.max_records or 200])),
            file=sys.stderr,
        )
    else:
        key_labels = np.ones(len(counts), dtype=int)
        order, labels, tree = list(range(n_records)), np.ones(n_records, dtype=int), None

    if args.scatter:
        with phase("project"):
            coordinates = embed_records(distances, args)
        with phase("scatter"):
            pairs = ((0, 1), (2, 3)) if args.embed == "pca" else ((0, 1),)
            draw_scatter(
                args.scatter, coordinates, counts, key_labels, palette, clustered, args, pairs
            )
    if args.no_tracks:
        report_timings()
        return

    with_dendrogram = clustered and args.dendrogram
    with_cluster_bar = clustered and args.cluster_bar
    cmap = plt.get_cmap(args.cmap)
    figure_height, row_height = track_figure_height(n_records, args)
    with phase("render"):
        raster = render_tracks(records, order, row_height, args, cmap)

    if not (with_dendrogram or with_cluster_bar or args.labels):
        # Nothing to draw around the tracks: the raster is the figure, and
        # writing it directly avoids building a canvas that, at tens of
        # thousands of rows, costs minutes and gigabytes inside matplotlib.
        with phase("write"):
            plt.imsave(args.output, raster)
        print("Saved %d records to %s" % (n_records, args.output), file=sys.stderr)
        report_timings()
        return
    figure = plt.figure(
        figsize=(args.width, figure_height),
        dpi=args.dpi,
        layout="constrained",  # Reserves room for the row labels, when they are drawn.
    )
    if with_dendrogram:
        grid = figure.add_gridspec(1, 2, width_ratios=[0.13, 1], wspace=0.01)
        ax_tree, ax = figure.add_subplot(grid[0, 0]), figure.add_subplot(grid[0, 1])
    else:
        ax_tree, ax = None, figure.add_subplot(1, 1, 1)

    # Row `position` of `order` occupies [position,position+1] vertically, and
    # the y axis points down, so that `order[0]` is at the top.
    ax.imshow(raster, extent=(0.0, 1.0, float(n_records), 0.0), aspect="auto", zorder=1)
    row_of = {leaf: position + 0.5 for position, leaf in enumerate(order)}
    if with_cluster_bar:  # Color bar of the flat clusters, left of the tracks.
        for leaf in order:
            ax.add_patch(
                Rectangle(
                    (-0.038, row_of[leaf] - 0.26),
                    0.018,
                    0.52,
                    facecolor=palette[(labels[leaf] - 1) % len(palette)],
                    edgecolor="none",
                    clip_on=False,
                )
            )

    ax.set_xlim(-0.045 if with_cluster_bar else 0.0, 1.0)
    ax.set_ylim(n_records, 0.0)
    ax.set_xticks(np.linspace(0, 1, 6))
    if args.labels:
        ax.set_yticks([row_of[leaf] for leaf in order])
        ax.set_yticklabels(
            [
                "%s%s  %s:%s  [%s:%d-%d]"
                % (
                    ("c%d  " % labels[leaf]) if clustered else "",
                    records[leaf]["file"],
                    records[leaf]["chrom"],
                    records[leaf]["pos"],
                    records[leaf]["coords"][0],
                    records[leaf]["coords"][1],
                    records[leaf]["coords"][2],
                )
                for leaf in order
            ],
            fontsize=5,
        )
        if with_dendrogram:  # Labels move right, to leave room for the tree.
            ax.yaxis.tick_right()
    else:
        ax.set_yticks([])
    ax.set_xlabel("Position in ALT (normalized)")
    ax.tick_params(axis="x", labelsize=7)
    ax.tick_params(axis="y", length=0, pad=2)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    if with_dendrogram:
        # The tree is over the keys, so a leaf sits at the middle of the block
        # of rows of the records that share that key.
        key_row, position = {}, 0
        for key in key_order:
            key_row[key] = position + len(members[key]) / 2.0
            position += len(members[key])
        maximum = draw_dendrogram(ax_tree, tree, key_row, key_labels, palette)
        ax_tree.set_xlim(maximum * 1.05, 0)  # Root on the left, leaves on the right.
        ax_tree.set_ylim(n_records, 0.0)
        ax_tree.set_yticks([])
        ax_tree.tick_params(axis="x", labelsize=6)
        ax_tree.set_xlabel("Distance", fontsize=7)
        for side in ("top", "right", "left"):
            ax_tree.spines[side].set_visible(False)

    with phase("write"):
        figure.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    plt.close(figure)
    print("Saved %d records to %s" % (n_records, args.output), file=sys.stderr)
    report_timings()


if __name__ == "__main__":
    main()
