#!/usr/bin/env python3
"""
Thread a long read onto a compact de Bruijn graph (GFA with sequences, e.g. built
by Bifrost) and emit a Bandage/BandageNG colour CSV highlighting the path.

Method: index every k-mer of every unitig (both orientations), then walk the read
k-mer by k-mer and record which unitig (and in which orientation) each k-mer
falls into.  Consecutive k-mers that stay inside the same unitig at consecutive
offsets are collapsed into a single visit.  Each visit becomes one step of the
path; the steps are coloured with a blue -> red gradient so that the direction of
travel is visible in Bandage.

Usage:
    ./thread_read_to_gfa.py graph.gfa reads.fa[.gz] -o colours.csv
    ./thread_read_to_gfa.py graph.gfa reads.fq.gz --read-id m64012_1234/567/ccs -o colours.csv
    ./thread_read_to_gfa.py graph.gfa --sequence ACGT... -o colours.csv

Then in Bandage: File -> Load CSV data, and switch the node colouring mode to
"Colour by CSV / custom colours".

Positional/optional arguments are described by --help.  k is taken from the GFA
header (Bifrost writes KL:Z:<k>) unless -k is given.
"""

import argparse
import colorsys
import csv
import gzip
import sys

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
ACGT = set("ACGT")

UNVISITED_COLOUR = "#d9d9d9"  # light grey for nodes the read does not touch


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def smart_open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_gfa(path):
    """Return (ordered list of segment names, {name: sequence}, k or None)."""
    names = []
    seqs = {}
    k = None
    with smart_open(path) as f:
        for line in f:
            if not line or line[0] not in "HS":
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "H":
                for tag in fields[1:]:
                    # Bifrost: KL:Z:63 ; some writers use KL:i:63
                    if tag.startswith("KL:Z:") or tag.startswith("KL:i:"):
                        k = int(tag[5:])
            elif fields[0] == "S":
                name, seq = fields[1], fields[2]
                if seq == "*":
                    sys.exit("ERROR: segment %s has no sequence ('*'); this "
                             "script needs a GFA with sequences." % name)
                names.append(name)
                seqs[name] = seq.upper()
    if not seqs:
        sys.exit("ERROR: no S lines found in %s" % path)
    return names, seqs, k


def iter_reads(path):
    """Yield (name, sequence) from a FASTA or FASTQ file (optionally gzipped)."""
    with smart_open(path) as f:
        first = f.readline()
        if not first:
            return
        if first[0] == ">":
            name = first[1:].strip().split()[0]
            chunks = []
            for line in f:
                if line[0] == ">":
                    yield name, "".join(chunks).upper()
                    name = line[1:].strip().split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip())
            yield name, "".join(chunks).upper()
        elif first[0] == "@":
            name = first[1:].strip().split()[0]
            while True:
                seq = f.readline().strip().upper()
                f.readline()  # '+'
                f.readline()  # qualities
                yield name, seq
                header = f.readline()
                if not header:
                    return
                name = header[1:].strip().split()[0]
        else:
            sys.exit("ERROR: %s is neither FASTA nor FASTQ" % path)


# ---------------------------------------------------------------------------
# k-mer index and threading
# ---------------------------------------------------------------------------

def build_kmer_index(seqs, k):
    """{kmer: (segment, offset, orientation)} for both orientations.

    orientation '+' means the k-mer is read off the segment left-to-right at
    `offset`, '-' means it is the reverse complement of that k-mer.
    """
    index = {}
    palindromic = 0
    for name, seq in seqs.items():
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            if not set(kmer) <= ACGT:
                continue
            rc = revcomp(kmer)
            if kmer == rc:
                palindromic += 1
            index[kmer] = (name, i, "+")
            index.setdefault(rc, (name, i, "-"))
    return index, palindromic


def thread_read(read_seq, index, k):
    """Walk the read and return (visits, n_kmers, n_matched).

    A visit is a dict describing one maximal run of read k-mers that stay in one
    segment in one orientation at consecutive offsets:
        segment, orient, read_start, read_end (inclusive k-mer indices),
        seg_start, seg_end (inclusive k-mer offsets in the segment).
    """
    visits = []
    n_kmers = max(0, len(read_seq) - k + 1)
    n_matched = 0
    current = None

    for i in range(n_kmers):
        hit = index.get(read_seq[i:i + k])
        if hit is None:
            current = None  # a gap always breaks the run
            continue
        n_matched += 1
        name, offset, orient = hit
        step = 1 if orient == "+" else -1
        if (current is not None
                and current["segment"] == name
                and current["orient"] == orient
                and current["read_end"] == i - 1
                and current["_next_offset"] == offset):
            current["read_end"] = i
            current["seg_end"] = offset
            current["_next_offset"] = offset + step
        else:
            current = {
                "segment": name,
                "orient": orient,
                "read_start": i,
                "read_end": i,
                "seg_start": offset,
                "seg_end": offset,
                "_next_offset": offset + step,
            }
            visits.append(current)

    for v in visits:
        v.pop("_next_offset")
    return visits, n_kmers, n_matched


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def gradient_colour(fraction):
    """Blue (start of the path) -> red (end of the path), as #RRGGBB."""
    hue = 0.66 * (1.0 - fraction)  # 0.66 = blue, 0.0 = red
    r, g, b = colorsys.hsv_to_rgb(hue, 0.85, 0.95)
    return "#%02x%02x%02x" % (int(r * 255), int(g * 255), int(b * 255))


def write_csv(out_path, seg_names, visits, include_unvisited, read_name):
    n = len(visits)
    per_segment = {}  # segment -> (colour of first visit, [step indices])
    for step, v in enumerate(visits):
        colour = gradient_colour(step / (n - 1) if n > 1 else 0.0)
        entry = per_segment.setdefault(v["segment"], [colour, []])
        entry[1].append(step + 1)

    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Name", "Colour", "Label", "PathSteps", "Read"])
        for name in seg_names:
            if name in per_segment:
                colour, steps = per_segment[name]
                label = "#" + ",".join(str(s) for s in steps)
                w.writerow([name, colour, label, len(steps), read_name])
            elif include_unvisited:
                w.writerow([name, UNVISITED_COLOUR, "", 0, ""])
    return per_segment


def main():
    ap = argparse.ArgumentParser(
        description="Thread a read onto a GFA de Bruijn graph and emit a "
                    "Bandage colour CSV.")
    ap.add_argument("gfa", help="GFA file with sequences (may be .gz)")
    ap.add_argument("reads", nargs="?",
                    help="FASTA/FASTQ with the read(s) (may be .gz)")
    ap.add_argument("--read-id", default=None,
                    help="ID of the read to thread (default: the first read)")
    ap.add_argument("--sequence", default=None,
                    help="thread this literal sequence instead of reading a file")
    ap.add_argument("-k", type=int, default=None,
                    help="k-mer size (default: from the GFA header KL tag)")
    ap.add_argument("-o", "--out", default="bandage_colours.csv",
                    help="output CSV (default: bandage_colours.csv)")
    ap.add_argument("--no-unvisited", action="store_true",
                    help="do not emit grey rows for segments off the path")
    ap.add_argument("--path-out", default=None,
                    help="also write the path (e.g. 12+,7-,...) to this file")
    args = ap.parse_args()

    if args.sequence is None and args.reads is None:
        ap.error("give a reads file or --sequence")

    seg_names, seqs, k_header = parse_gfa(args.gfa)
    k = args.k or k_header
    if k is None:
        ap.error("k not found in the GFA header; pass -k")
    print("Graph: %d segments, k=%d" % (len(seg_names), k), file=sys.stderr)

    if args.sequence is not None:
        read_name, read_seq = "sequence", args.sequence.upper()
    else:
        read_seq = None
        for name, seq in iter_reads(args.reads):
            if args.read_id is None or name == args.read_id:
                read_name, read_seq = name, seq
                break
        if read_seq is None:
            sys.exit("ERROR: read %s not found in %s" % (args.read_id, args.reads))
    print("Read:  %s (%d bp)" % (read_name, len(read_seq)), file=sys.stderr)
    if len(read_seq) < k:
        sys.exit("ERROR: read is shorter than k")

    index, palindromic = build_kmer_index(seqs, k)
    if palindromic:
        print("Warning: %d palindromic k-mers (orientation is ambiguous there)"
              % palindromic, file=sys.stderr)

    visits, n_kmers, n_matched = thread_read(read_seq, index, k)
    print("k-mers: %d in the read, %d found in the graph (%.2f%%)"
          % (n_kmers, n_matched, 100.0 * n_matched / n_kmers if n_kmers else 0.0),
          file=sys.stderr)
    print("Path:   %d steps over %d distinct segments"
          % (len(visits), len({v["segment"] for v in visits})), file=sys.stderr)

    path_str = ",".join(v["segment"] + v["orient"] for v in visits)
    if args.path_out:
        with open(args.path_out, "w") as fh:
            fh.write(">%s\n%s\n" % (read_name, path_str))
    else:
        print(path_str)

    # Report the gaps: places where the read leaves the graph.
    prev_end = -1
    for v in visits:
        if v["read_start"] > prev_end + 1:
            print("  gap: read k-mers %d-%d unmatched"
                  % (prev_end + 1, v["read_start"] - 1), file=sys.stderr)
        prev_end = v["read_end"]
    if n_kmers and prev_end < n_kmers - 1:
        print("  gap: read k-mers %d-%d unmatched" % (prev_end + 1, n_kmers - 1),
              file=sys.stderr)

    write_csv(args.out, seg_names, visits, not args.no_unvisited, read_name)
    print("Wrote %s" % args.out, file=sys.stderr)


if __name__ == "__main__":
    main()
