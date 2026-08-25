#!/usr/bin/env python3
"""Count BNDs per sample in a multi-sample VCF, split by the kind of contig pair
they join, and save the categories as a jittered scatterplot.

The name of every sample with at least `min_intra` intra-chromosomal BNDs (of any
of the three intra-chromosomal categories below) or at least `min_inter`
inter-chromosomal ones is printed to stdout, one per line, in the order the
samples appear in the VCF header; the two thresholds are ORed, so crossing
either is enough.

One point per sample per category: its y is the number of BND records that are
non-ref (het or hom-alt) in that sample, its x is the category position plus a
small random offset, so that samples with equal counts stay distinguishable.

A second panel, to the right of the first one and drawn the same way, breaks the
same records down by standard chromosome instead: one column per chr1...chr22,
chrX, chrY, chrM, and a record counts for each of the (one or two) standard
contigs it touches. The two panels have independent y axes.

A contig is *standard* when it is one of chr1...chr22, chrX, chrY, chrM (names
without the `chr` prefix and the `MT` alias are accepted too); everything else
-- `chr*_*_random`, `chrUn_*`, `chr*_alt`, `chrEBV`, decoys -- is non-standard.
The categories **overlap**: each record is counted in every one it matches, so
the columns do not add up to the total number of BNDs.

  intra DEL-like                both contigs the same, deletion orientation
  intra DUP-like                both contigs the same, duplication orientation
  intra INV-like                both contigs the same, any other orientation
  inter-chromosomal             the two contigs differ
  chrM + standard               one contig chrM, the other a different standard
  chrEBV + standard             one contig chrEBV, the other standard
  standard + non-standard       exactly one contig non-standard
  non-standard + non-standard   both contigs non-standard, same one included

A chrEBV-to-chr1 breakend, for instance, lands in inter-chromosomal, in
chrEBV + standard and in standard + non-standard at once.

An intra-chromosomal record at POS=X with mate position Y is DEL-like when its
ALT is `X[chr:Y[` with Y>X or `]chr:Y]X` with Y<X, and DUP-like when it is
`X[chr:Y[` with Y<X or `]chr:Y]X` with Y>X. Everything else -- the two remaining
orientations `X]chr:Y]` and `[chr:Y[X`, Y==X, and records whose mate position is
not known because the mate came from INFO/CHR2 rather than from a bracketed ALT
-- is INV-like.

After the sample names come the cohort-wide record counts -- intra-chromosomal,
inter-chromosomal, and the three intra-chromosomal sub-categories -- counted per
record rather than per sample, so a record counts once however many samples carry
it, and then the widest three records of each of those three sub-categories,
widest first, where the span of a record is |Y-X|. Every line of that report is either a '#' comment -- naming the
category, or giving the span and carrier count of the record below it -- or the
record itself: its nine non-sample columns verbatim, with the whole sample block
replaced by a comma-separated list of the samples that have a non-ref genotype
for it. Only records with at least one such sample are eligible, and only those
with a known mate position: the INV-like records that got there by having none
have no span, so they are counted but never reported.

The VCF is scanned once, top to bottom. It is read as bytes and the per-record
sample block is never split into Python strings: genotypes are scanned with
numpy at C speed, which keeps cohort-scale VCFs (>10k samples) practical.
"""

import sys
import gzip
import heapq
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")          # write-to-file only; no display needed
import matplotlib.pyplot as plt

TAB = 9
COLON = 58
DIGIT_1 = 49
DIGIT_9 = 57

DEL_LIKE = 0
DUP_LIKE = 1
INV_LIKE = 2
INTER = 3
M_STD = 4
EBV_STD = 5
STD_NONSTD = 6
NONSTD_NONSTD = 7
N_CATEGORIES = 8

# The three ways an intra-chromosomal record can be classified; a record lands in
# exactly one of them, so summing these rows gives the intra-chromosomal total.
INTRA_CATEGORIES = (DEL_LIKE, DUP_LIKE, INV_LIKE)

# Single-line names, used for the TSV header.
CATEGORIES = ["intra DEL-like", "intra DUP-like", "intra INV-like",
              "inter-chromosomal", "chrM + standard", "chrEBV + standard",
              "standard + non-standard", "non-standard + non-standard"]
# Same order, wrapped so that eight ticks fit under the axis.
CATEGORY_LABELS = ["intra\nDEL-like", "intra\nDUP-like", "intra\nINV-like",
                   "inter-\nchromosomal", "chrM +\nstandard",
                   "chrEBV +\nstandard", "standard +\nnon-standard",
                   "non-standard +\nnon-standard"]
# The three intra-chromosomal columns keep the blue of the single intra column
# they replace, as three shades, so that they still read as one family.
COLORS = ["#1f4e79", "#2f6db3", "#7fa9d4", "#c2650f", "#3f8f5b", "#a04ba0",
          "#b3312f", "#7a6a55"]

# How many widest-span records per intra-chromosomal category to print.
N_LONGEST = 3

# Which side of the inserted sequence the bracket pair sits on: `t[chr:Y[` and
# `t]chr:Y]` are BRACKET_AFTER, `[chr:Y[t` and `]chr:Y]t` are BRACKET_BEFORE.
BRACKET_AFTER = 0
BRACKET_BEFORE = 1

MITO = b"chrM"
EBV = b"chrEBV"
# Column order of the per-chromosome panel; also the definition of "standard".
CHROMOSOMES = ([b"chr%d" % i for i in range(1, 23)]
               + [b"chrX", b"chrY", MITO])
CHROMOSOME_INDEX = {name: i for i, name in enumerate(CHROMOSOMES)}
N_CHROMOSOMES = len(CHROMOSOMES)
# The `chr` prefix is dropped: 25 ticks only fit if the labels are short.
CHROMOSOME_LABELS = [name.decode()[3:] for name in CHROMOSOMES]
# Alternated column by column, which is all that is needed to tell the narrow
# neighbouring clouds apart.
CHROMOSOME_COLORS = ["#2f6db3", "#79a8d8"]
STANDARD = frozenset(CHROMOSOMES)


def canonical(contig):
    """Return `contig` in `chr`-prefixed form, with `chrMT` folded onto chrM.

    Callsets built against the no-prefix flavour of GRCh38 name the contigs
    `1`, `X`, `MT`; normalising here lets one contig table serve both.
    """
    name = contig if contig.startswith(b"chr") else b"chr" + contig
    return MITO if name == b"chrMT" else name


def intra_category(pos, mate_pos, orientation):
    """Classify an intra-chromosomal BND as DEL-, DUP- or INV-like.

    Of the four breakend orientations only two describe a plain deletion or
    tandem duplication, and which of the two they are depends on whether the
    mate lies to the right or to the left of POS:

      X[chr:Y[   Y>X deletion-like,    Y<X duplication-like
      ]chr:Y]X   Y<X deletion-like,    Y>X duplication-like

    The other two, X]chr:Y] and [chr:Y[X, join two same-strand-broken pieces and
    are inversion-like whatever Y is. Y==X, and a mate position that is not known
    at all (`orientation` or `mate_pos` None, i.e. the mate came from INFO/CHR2),
    also land in INV-like: it is the catch-all category.
    """
    if pos is None or mate_pos is None or orientation is None or mate_pos == pos:
        return INV_LIKE
    side, bracket = orientation
    if side == BRACKET_AFTER and bracket == b"[":
        return DEL_LIKE if mate_pos > pos else DUP_LIKE
    if side == BRACKET_BEFORE and bracket == b"]":
        return DEL_LIKE if mate_pos < pos else DUP_LIKE
    return INV_LIKE


def categories_of(chrom, pos, mate, mate_pos, orientation):
    """Return every category index a BND joining `chrom` and `mate` falls into.

    The categories overlap, so this is a list, never a single index. The
    intra-chromosomal trio and inter-chromosomal together always contribute
    exactly one entry; how standard the two contigs are contributes zero, one or
    two more.
    """
    first = canonical(chrom)
    second = canonical(mate)
    first_standard = first in STANDARD
    second_standard = second in STANDARD

    found = [intra_category(pos, mate_pos, orientation)
             if first == second else INTER]
    if not first_standard and not second_standard:
        found.append(NONSTD_NONSTD)
    elif not first_standard or not second_standard:
        found.append(STD_NONSTD)
        if (first if not first_standard else second) == EBV:
            found.append(EBV_STD)
    elif (first == MITO) != (second == MITO):
        # Both standard and exactly one is chrM, so the mate is "another"
        # standard chromosome; a chrM-to-chrM breakend is intra-chromosomal.
        found.append(M_STD)
    return found


def chromosomes_of(chrom, mate):
    """Return the column indices of the standard contigs a BND touches.

    Zero, one or two entries: non-standard contigs have no column, and an
    intra-chromosomal breakend touches its chromosome once, not twice.
    """
    found = []
    for contig in {canonical(chrom), canonical(mate)}:
        column = CHROMOSOME_INDEX.get(contig)
        if column is not None:
            found.append(column)
    return found


def open_vcf(path):
    """Open a plain or bgzipped/gzipped VCF in binary mode; '-' means stdin."""
    if path == "-":
        return sys.stdin.buffer
    if path.endswith(".gz") or path.endswith(".bgz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def mate_of(alt, info):
    """Return (contig, position, orientation) of a BND's mate.

    The breakend ALT syntax puts the mate locus between a pair of brackets --
    't[chr2:1000[', ']chr2:1000]t' and the other two orientations -- so the
    contig is what precedes the last ':' inside the brackets and the position is
    what follows it. `orientation` is the (side, bracket) pair the classification
    in `intra_category` needs: which side of the inserted sequence the brackets
    sit on, and which of the two bracket characters they are.

    Symbolic <BND> records without brackets fall back to INFO/CHR2 for the
    contig, and get a None position and orientation; single breakends ('.t',
    't.') have no mate at all and yield (None, None, None).
    """
    first = alt.split(b",")[0]
    open_bracket = -1
    close_bracket = -1
    for bracket in (b"[", b"]"):
        position = first.find(bracket)
        if position >= 0 and (open_bracket < 0 or position < open_bracket):
            open_bracket = position
            close_bracket = first.find(bracket, position + 1)
    if open_bracket >= 0 and close_bracket > open_bracket:
        locus = first[open_bracket + 1:close_bracket]
        cut = locus.rfind(b":")
        if cut > 0:
            side = BRACKET_BEFORE if open_bracket == 0 else BRACKET_AFTER
            orientation = (side, first[open_bracket:open_bracket + 1])
            try:
                mate_pos = int(locus[cut + 1:])
            except ValueError:      # not a plain integer: treat as unknown
                mate_pos = None
            return locus[:cut], mate_pos, orientation

    start = info.find(b"CHR2=")
    while start >= 0:
        if start == 0 or info[start - 1] == ord(";"):   # a real INFO key
            end = info.find(b";", start)
            return info[start + 5:end if end >= 0 else len(info)], None, None
        start = info.find(b"CHR2=", start + 1)
    return None, None, None


def alt_mask(sample_block, n_samples):
    """Return a boolean array: True where that sample's GT has a non-ref allele.

    `sample_block` is the raw bytes of the record from the first sample column
    to the end of the line (trailing newline included, harmlessly). A sample is
    ALT iff its GT sub-field -- everything before the first ':' of the column,
    per the VCF spec, which requires GT to come first when present -- contains a
    digit 1-9. That covers '0/1', '1|1', '1/2', multi-digit allele indices and
    haploid calls, while '0/0', './.' and '.' stay false.

    The whole block is scanned with vector ops instead of one Python call per
    sample, which is what makes the 12k-sample case tractable.
    """
    array = np.frombuffer(sample_block, dtype=np.uint8)
    is_tab = array == TAB
    starts = np.empty(int(is_tab.sum()) + 1, dtype=np.int64)
    starts[0] = 0
    np.add(np.flatnonzero(is_tab), 1, out=starts[1:])
    if len(starts) != n_samples:
        raise RuntimeError("Record has %d sample columns, header declares %d."
                           % (len(starts), n_samples))

    # A character belongs to its column's GT sub-field iff no ':' has appeared
    # since the column started, i.e. the running colon count still equals the
    # count at the column's first character.
    colons = np.cumsum(array == COLON, dtype=np.int32)
    at_start = np.zeros(len(starts), dtype=np.int32)
    at_start[1:] = colons[starts[1:] - 1]
    column_of = np.cumsum(is_tab, dtype=np.int32)
    in_gt = colons == at_start[column_of]

    non_ref = (array >= DIGIT_1) & (array <= DIGIT_9) & in_gt
    return np.bitwise_or.reduceat(non_ref, starts)


def keep_longest(heap, span, number, columns, carriers):
    """Add a record to a bounded top-N-by-span min-heap, dropping the smallest.

    `heap` never grows past N_LONGEST, so the whole scan holds a handful of
    records rather than one entry per record. Records with no ALT sample are not
    candidates, and `number` -- unique per record -- breaks span ties so that the
    heap never has to order two numpy arrays.
    """
    if not len(carriers):
        return
    entry = (span, number, columns, carriers)
    if len(heap) < N_LONGEST:
        heapq.heappush(heap, entry)
    elif span > heap[0][0]:
        heapq.heapreplace(heap, entry)


def build_counts(vcf_path, progress_every=0):
    """Stream the VCF once, returning (counts, chrom_counts, totals, longest,
    samples, n_skipped).

    `counts` is an N_CATEGORIES x n_samples int64 array, one row per entry of
    CATEGORIES; a record contributes to every category it matches, so the rows
    overlap. `chrom_counts` is an N_CHROMOSOMES x n_samples int64 array in
    CHROMOSOMES order, holding the same records broken down by the standard
    contigs they touch. `n_skipped` counts records whose mate contig could not
    be determined (single breakends, or symbolic ALTs without INFO/CHR2); those
    contribute to neither array.

    `totals` is an N_CATEGORIES int64 array counting the same categories over
    *records* instead of over sample genotypes: one increment per record per
    category it matches, whether or not any sample carries it. It is what the
    cohort has, where `counts` is what each sample has.

    `longest` is a list indexed by intra-chromosomal category, each holding the
    N_LONGEST widest-span records of that category as a min-heap of
    (span, record number, the nine non-sample columns, indices of the ALT
    samples) tuples. Only records with at least one ALT sample are candidates,
    for the same reason the counts are per-sample and non-ref: a record no
    sample carries says nothing about the cohort.
    """
    samples = None
    counts = None
    chrom_counts = None
    totals = np.zeros(N_CATEGORIES, dtype=np.int64)
    longest = [[] for _ in INTRA_CATEGORIES]
    n_records = 0
    n_skipped = 0

    with open_vcf(vcf_path) as vcf:
        for line in vcf:
            if line.startswith(b"#"):
                if line.startswith(b"##"):
                    continue
                samples = line.rstrip().decode().split("\t")[9:]
                counts = np.zeros((N_CATEGORIES, len(samples)), dtype=np.int64)
                chrom_counts = np.zeros((N_CHROMOSOMES, len(samples)),
                                        dtype=np.int64)
                continue
            if counts is None:
                raise RuntimeError("Malformed VCF: no #CHROM header line.")

            # Split off only the 9 leading columns; the sample block stays one
            # bytes object rather than becoming one Python string per sample.
            columns = line.split(b"\t", 9)
            if len(columns) < 10:
                continue
            mate, mate_pos, orientation = mate_of(columns[4], columns[7])
            if mate is None:
                n_skipped += 1
                continue
            try:
                pos = int(columns[1])
            except ValueError:      # unparsable POS: no DEL/DUP/INV call
                pos = None
            # One genotype scan per record, added to each category it matches.
            mask = alt_mask(columns[9], len(samples))
            found = categories_of(columns[0], pos, mate, mate_pos, orientation)
            for category in found:
                counts[category] += mask
                totals[category] += 1
            for column in chromosomes_of(columns[0], mate):
                chrom_counts[column] += mask
            # found[0] is the intra-or-inter entry, so this is the record's
            # intra-chromosomal category when it has one.
            if found[0] in INTRA_CATEGORIES and pos is not None \
                    and mate_pos is not None:
                keep_longest(longest[found[0]], abs(mate_pos - pos), n_records,
                             columns[:9], np.flatnonzero(mask))

            n_records += 1
            if progress_every and n_records % progress_every == 0:
                sys.stderr.write("%d records counted\n" % n_records)
                sys.stderr.flush()

    if counts is None:
        raise RuntimeError("Malformed VCF: no #CHROM header line.")
    return counts, chrom_counts, totals, longest, samples, n_skipped


def draw_clouds(axes, matrix, colors, generator, size, alpha):
    """Scatter one jittered cloud per row of `matrix`, and return the medians.

    Column `i` of the panel holds row `i` of the matrix, drawn at x=i plus a
    random offset per sample; a short horizontal rule marks the median, the one
    summary worth reading off a cloud directly.
    """
    n_samples = matrix.shape[1]
    medians = []
    for column in range(matrix.shape[0]):
        values = matrix[column]
        jitter = generator.uniform(-0.18, 0.18, size=n_samples)
        axes.scatter(column + jitter, values, s=size, alpha=alpha,
                     color=colors[column % len(colors)], linewidths=0.0,
                     zorder=2)
        median = float(np.median(values)) if n_samples else 0.0
        axes.hlines(median, column - 0.3, column + 0.3, color="#333333",
                    linewidth=2.0, zorder=3)
        medians.append(median)
    return medians


def style_panel(axes, n_columns, log_scale):
    """Apply the shared look: y grid behind the points, no top/right spines."""
    axes.set_xlim(-0.6, n_columns - 0.4)
    if log_scale:
        axes.set_yscale("log")
    else:
        axes.set_ylim(bottom=0)
    axes.grid(axis="y", color="#dddddd", linewidth=0.8, zorder=0)
    axes.set_axisbelow(True)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)


def plot_counts(counts, chrom_counts, samples, image_path, dpi, seed,
                log_scale):
    """Draw the category panel and the per-chromosome one, write `image_path`.

    Both panels are jittered point clouds of the same quantity -- non-ref BND
    records per sample -- binned differently on x. Their y axes are independent:
    inter-chromosomal alone outweighs any single chromosome by enough that a
    shared axis would flatten the right panel into a band.
    The format is taken from the file extension (.png, .pdf, .svg, ...).
    """
    n_samples = len(samples)
    generator = np.random.default_rng(seed)
    figure, (left, right) = plt.subplots(
        1, 2, figsize=(19.0, 6.0), gridspec_kw={"width_ratios": [4.0, 4.5]})

    # Denser clouds need smaller, more transparent points to stay readable.
    if n_samples <= 100:
        size, alpha = 42.0, 0.85
    elif n_samples <= 2000:
        size, alpha = 16.0, 0.5
    else:
        size, alpha = 6.0, 0.25

    medians = draw_clouds(left, counts, COLORS, generator, size, alpha)
    # With eight clouds side by side there is no room for a label next to the
    # median line, so it goes under the category name instead.
    left.set_xticks(list(range(N_CATEGORIES)))
    left.set_xticklabels(["%s\nmedian %.0f" % (label, median)
                          for label, median in zip(CATEGORY_LABELS, medians)],
                         fontsize=8)
    left.set_xlabel("BND category (categories overlap; a record can be in several)")
    left.set_ylabel("number of non-ref BND records per sample")
    style_panel(left, N_CATEGORIES, log_scale)

    # 25 columns leave room for the chromosome name and nothing else, so the
    # medians here are only the drawn rules.
    draw_clouds(right, chrom_counts, CHROMOSOME_COLORS, generator, size, alpha)
    right.set_xticks(list(range(N_CHROMOSOMES)))
    right.set_xticklabels(CHROMOSOME_LABELS, fontsize=8)
    right.set_xlabel("chromosome (a record counts once for each standard "
                     "contig it touches)")
    right.set_ylabel("number of non-ref BND records per sample")
    style_panel(right, N_CHROMOSOMES, log_scale)

    figure.suptitle("BND calls per sample (n=%d)" % n_samples)
    figure.tight_layout()
    figure.savefig(image_path, dpi=dpi)
    plt.close(figure)


def outlier_samples(counts, samples, min_intra, min_inter):
    """Return the samples over either threshold, in VCF header order.

    The thresholds are ORed and inclusive: a sample is returned when its
    intra-chromosomal count -- the DEL-, DUP- and INV-like columns together,
    since every intra-chromosomal record is in exactly one of them -- is >=
    `min_intra`, or its inter-chromosomal one is >= `min_inter`.
    """
    intra = counts[list(INTRA_CATEGORIES)].sum(axis=0)
    selected = (intra >= min_intra) | (counts[INTER] >= min_inter)
    return [sample for sample, keep in zip(samples, selected) if keep]


def report_totals(totals):
    """Print the cohort-wide record count of the intra and inter categories.

    These are record counts, so a record counts once however many samples carry
    it -- unlike the plot and the TSV, which count per sample and only non-ref
    genotypes. The intra-chromosomal total is the sum of the three sub-categories
    below it, since every intra-chromosomal record is in exactly one of them.
    """
    intra = int(totals[list(INTRA_CATEGORIES)].sum())
    print("# intra-chromosomal BNDs\t%d" % intra)
    print("# inter-chromosomal BNDs\t%d" % int(totals[INTER]))
    for category in INTRA_CATEGORIES:
        print("# %s BNDs\t%d" % (CATEGORIES[category], int(totals[category])))


def report_longest(longest, samples):
    """Print the widest-span records of each intra-chromosomal category.

    One block per category, widest first: a '#' comment naming the category, then
    a '#' comment with the span and carrier count of each record followed by the
    record itself -- its nine non-sample columns verbatim, with the whole sample
    block replaced by a comma-separated list of the samples whose genotype is
    non-ref. Keeping the span out of the record line leaves that line a plain
    tab-separated VCF prefix.
    """
    for category in INTRA_CATEGORIES:
        # Equal spans fall back to VCF order, so the report is deterministic.
        records = sorted(longest[category],
                         key=lambda entry: (-entry[0], entry[1]))
        if not records:
            print("# %s: no records with a span and an ALT sample"
                  % CATEGORIES[category])
            continue
        print("# %s: %d widest-span record%s, widest first"
              % (CATEGORIES[category], len(records),
                 "" if len(records) == 1 else "s"))
        for span, _, columns, carriers in records:
            print("# span %d, %d ALT sample%s"
                  % (span, len(carriers), "" if len(carriers) == 1 else "s"))
            print("\t".join([column.decode() for column in columns]
                            + [",".join(samples[i] for i in carriers)]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="multi-sample BND-only VCF (plain or gzipped; - for stdin)")
    parser.add_argument("image", help="output plot file; the extension picks "
                                      "the format (.png, .pdf, .svg, ...)")
    parser.add_argument("min_intra", type=int, metavar="X",
                        help="print a sample if it has >=X intra-chromosomal BNDs")
    parser.add_argument("min_inter", type=int, metavar="Y",
                        help="print a sample if it has >=Y inter-chromosomal BNDs")
    parser.add_argument("--tsv", help="also write the per-sample counts to this TSV file")
    parser.add_argument("--log", action="store_true",
                        help="use a logarithmic y axis")
    parser.add_argument("--dpi", type=int, default=150,
                        help="resolution of the output image (default: 150)")
    parser.add_argument("--seed", type=int, default=42,
                        help="seed of the x jitter, for reproducible plots (default: 42)")
    parser.add_argument("--progress", type=int, default=0, metavar="N",
                        help="report to stderr every N counted records")
    args = parser.parse_args()

    counts, chrom_counts, totals, longest, samples, n_skipped = build_counts(
        args.vcf, args.progress)
    if n_skipped:
        sys.stderr.write("Warning: %d records skipped (no mate contig in ALT "
                         "or INFO/CHR2).\n" % n_skipped)
    for sample in outlier_samples(counts, samples, args.min_intra,
                                  args.min_inter):
        print(sample)
    report_totals(totals)
    report_longest(longest, samples)
    if args.tsv:
        with open(args.tsv, "wt") as out:
            out.write("sample\t" + "\t".join(CATEGORIES) + "\t"
                      + "\t".join(name.decode() for name in CHROMOSOMES) + "\n")
            for j, sample in enumerate(samples):
                out.write(sample + "\t"
                          + "\t".join(str(counts[c, j])
                                      for c in range(N_CATEGORIES)) + "\t"
                          + "\t".join(str(chrom_counts[c, j])
                                      for c in range(N_CHROMOSOMES)) + "\n")
    plot_counts(counts, chrom_counts, samples, args.image, args.dpi, args.seed,
                args.log)


if __name__ == "__main__":
    main()
