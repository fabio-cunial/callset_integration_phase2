#!/usr/bin/env python3
"""
Strip all sequence information from a GFA file, replacing every segment sequence
with a single artificial character.  Topology (links, containments, paths,
walks, jumps) and all optional tags are preserved verbatim, so the output is a
much smaller graph that can still be loaded by Bandage / gfatools / vg and
inspected structurally.

Both GFA1 (S <name> <seq> [tags]) and GFA2 (S <sid> <slen> <seq> [tags]) are
handled; the flavour is detected per S line from the shape of the record.

By default the original sequence length is kept in an LN:i tag (added if
missing, overwritten if already present) so that no length information is lost.
Note that after stripping, the segment lengths implied by the sequences no
longer agree with the CIGARs on L/C/E lines; pass --clear-overlaps to replace
those with "*" if a downstream tool complains.

Usage:
    ./strip_gfa_sequences.py graph.gfa stripped.gfa
    ./strip_gfa_sequences.py graph.gfa.gz stripped.gfa.gz -c N
    ./strip_gfa_sequences.py - - < graph.gfa > stripped.gfa

Input/output may be "-" for stdin/stdout; ".gz" paths are (de)compressed.
"""

import argparse
import gzip
import sys


def smart_open_read(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def smart_open_write(path):
    if path == "-":
        return sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "wt")


def is_tag(field):
    """True if the field looks like an optional tag TAG:T:VALUE."""
    return (
        len(field) >= 5
        and field[2] == ":"
        and field[4] == ":"
        and field[0].isalpha()
        and field[1].isalnum()
    )


def set_ln_tag(tags, length):
    """Return tags with LN:i:<length> overwritten in place, or appended if absent."""
    out = list(tags)
    for i, t in enumerate(out):
        if t.startswith("LN:i:"):
            out[i] = "LN:i:%d" % length
            return out
    out.append("LN:i:%d" % length)
    return out


def strip_segment(fields, char, keep_length):
    """Rewrite one S line (given as a list of fields) in place-ish; returns fields."""
    # GFA2:  S <sid> <slen> <sequence> [tags]     -> field 2 is an integer length
    # GFA1:  S <name> <sequence> [tags]
    gfa2 = len(fields) >= 4 and fields[2].isdigit() and not is_tag(fields[3])
    seq_idx = 3 if gfa2 else 2
    if len(fields) <= seq_idx:
        return fields  # malformed / truncated S line: leave alone

    seq = fields[seq_idx]
    tags = fields[seq_idx + 1:]

    if seq == "*":
        # No sequence to strip; still normalise the length if we know it.
        length = int(fields[2]) if gfa2 else None
    else:
        length = len(seq)
        fields[seq_idx] = char

    if gfa2:
        fields[2] = "1" if seq != "*" else fields[2]
    if keep_length and length is not None:
        tags = set_ln_tag(tags, length)
    fields[seq_idx + 1:] = tags
    return fields


def main():
    ap = argparse.ArgumentParser(
        description="Strip sequences from a GFA, one artificial character per node.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("input", help='input GFA ("-" for stdin, ".gz" allowed)')
    ap.add_argument("output", help='output GFA ("-" for stdout, ".gz" allowed)')
    ap.add_argument(
        "-c", "--char", default="N",
        help="the single artificial character to use for every node sequence",
    )
    ap.add_argument(
        "--no-length", action="store_true",
        help="do not add/update LN:i tags (existing ones are still passed through)",
    )
    ap.add_argument(
        "--clear-overlaps", action="store_true",
        help='replace the overlap/CIGAR field of L, C and E lines with "*"',
    )
    args = ap.parse_args()

    if len(args.char) != 1:
        ap.error("--char must be exactly one character, got %r" % args.char)

    n_segments = 0
    n_stripped_bp = 0
    fin = smart_open_read(args.input)
    fout = smart_open_write(args.output)
    try:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            rec = line[0]
            if rec == "S":
                fields = line.split("\t")
                before = len(line)
                fields = strip_segment(fields, args.char, not args.no_length)
                line = "\t".join(fields)
                n_segments += 1
                n_stripped_bp += max(0, before - len(line))
            elif args.clear_overlaps and rec in "LCE":
                fields = line.split("\t")
                # GFA1: L <from> <fo> <to> <to> <overlap>; C <c> <co> <ct> <to> <pos> <overlap>
                # GFA2: E <eid> <sid1> <sid2> <b1> <e1> <b2> <e2> <alignment>
                ovl_idx = {"L": 5, "C": 6, "E": 8}[rec]
                if len(fields) > ovl_idx and not is_tag(fields[ovl_idx]):
                    fields[ovl_idx] = "*"
                    line = "\t".join(fields)
            fout.write(line)
            fout.write("\n")
    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

    sys.stderr.write(
        "Stripped %d segments, ~%d characters of sequence removed.\n"
        % (n_segments, n_stripped_bp)
    )


if __name__ == "__main__":
    main()
