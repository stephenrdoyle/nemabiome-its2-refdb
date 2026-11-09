#!/usr/bin/env python3
"""The ITS2 reference unit: re-cut every backbone record at fixed anchors.

    bin/its2.py <seqs_nr.fasta> <records.fasta.gz> <out prefix>

Reads from the NC1/NC2 amplicon are 5.8S tail + ITS2 + ~47 bases of 28S.
A reference must cover the same stretch as the reads it is compared with,
or identities mean different things from record to record. The nemabiome
backbone's ITS2 (cut by markerDB's HMM) does not: where a source record
runs to NC2 the HMM stops before the 28S, where it ends sooner the HMM
keeps the 28S stretch, and the start wanders by up to 20 bases.

So the unit is fixed by two anchors, not by the HMM:

  start  the 3' end of NC1 (the 5.8S end; ITS2 begins 3 bases earlier,
         inside the primer, which reads lose with primer trimming anyway)
  end    the start of the 28S stretch before NC2 (S28 below)

and amphid scores each read only over the span its reference covers, so
the read's own 28S stretch counts neither for nor against. The 28S
stretch varies by 1-6 bases between nematodes; borrowing one onto a
record that lacks it would invent differences larger than the species
margin, so it is never added.

Each end is located, in order:
  primer    NC1 within 3 edits / S28 within 8 edits, anywhere near the
            backbone boundary
  loose     within 6 / 12 edits, in a window around the backbone
            boundary (taxa outside the Strongylida, whose primer site
            differs by 4-6 bases)
  record    the record itself starts / ends at the backbone boundary:
            the submitter trimmed there, so the end is unanchored
  hmm       flanking sequence present but no anchor found: the
            backbone's boundary is kept (start moved 3 bases in, to
            where NC1 would end)

Writes <prefix>.fasta (re-cut ITS2, forward strand, keyed by accession)
and <prefix>.tsv (accession, start/end method and edit cost, length,
the 28S stretch where the record runs to NC2, and NC1/NC2 edits where
the record holds the primer site - blank where it does not).
"""

from __future__ import annotations

import collections
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seqmatch import find_anchor_cost, revcomp  # noqa: E402

# Gasser et al. 1993 (Nucleic Acids Res 21:2525), as used by Avramenko
# et al. 2015 for the nemabiome. NC1 sits at the 3' end of 5.8S, NC2 in 28S.
NC1 = "ACGTCTGGTTCAGGGTTGTT"
NC2 = "TTAGTTTCTTTTCCTCCGCT"
# The head of the 28S stretch between ITS2 and NC2: the commonest of 20
# variants in 109 GenBank records spanning NC1..NC2, less its last 3 bases.
S28 = "AACCTGAGCTCAGTCGTGATTACCCGCCGAACTTAAGCATATCA"
TIGHT = {"NC1": 3, "S28": 8}
LOOSE = {"NC1": 6, "S28": 12}
WINDOW = 30         # loose search: bases either side of the backbone boundary
HMM_START_SHIFT = 3  # backbone ITS2 starts 3 bases before NC1 ends (2,544 of 3,820)


def read_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, cur, buf = {}, None, []
    for line in op(path, "rt"):
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur, buf = line[1:].split()[0].split(".")[0], []
        else:
            buf.append(line.strip().upper())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def locate(record, its):
    """The record on the backbone's strand, and the backbone's span in it."""
    for x in (record, revcomp(record)):
        p = x.find(its[:30])
        q = x.rfind(its[-30:])
        if p >= 0 and q >= p:
            return x, p, q + 30
    return None, -1, -1


def cut_start(x, p):
    # Long records run back through 5.8S and ITS1 into 18S: stay near.
    lo = max(0, p - 90)
    hit = find_anchor_cost(x[lo:p + WINDOW], NC1, TIGHT["NC1"])
    if hit:
        return lo + hit[1], "primer", hit[2]
    lo = max(0, p - WINDOW)
    hit = find_anchor_cost(x[lo:p + WINDOW], NC1, LOOSE["NC1"])
    if hit:
        return lo + hit[1], "loose", hit[2]
    if p == 0:
        return 0, "record", ""
    return p + HMM_START_SHIFT, "hmm", ""


def cut_end(x, q):
    # The backbone's end may sit before the 28S stretch or after it (where
    # the HMM swallowed it), so the window reaches a stretch's length past.
    lo = max(0, q - 60)
    hit = find_anchor_cost(x[lo:q + 60], S28, TIGHT["S28"])
    method = "primer"
    if not hit:
        lo = max(0, q - WINDOW - 47)
        hit = find_anchor_cost(x[lo:q + WINDOW], S28, LOOSE["S28"])
        method = "loose"
    if hit:
        return lo + hit[0], method, hit[2]
    if q == len(x):
        return q, "record", ""
    return q, "hmm", ""


def main() -> int:
    bb_path, rec_path, prefix = sys.argv[1:4]
    bb, rec = read_fasta(bb_path), read_fasta(rec_path)
    rows, tally = [], collections.Counter()
    nc2r = revcomp(NC2)
    with open(prefix + ".fasta", "w") as fa:
        for acc, its in bb.items():
            if acc not in rec:
                tally["no record"] += 1
                continue
            x, p, q = locate(rec[acc], its)
            if x is None:
                tally["backbone not in record"] += 1
                continue
            s, sm, sc = cut_start(x, p)
            e, em, ec = cut_end(x, q)
            if e - s > len(its) + 30:
                # Anchors can only trim what the HMM over-reached, never add
                # much: a cut far longer than the backbone ITS2 found an
                # anchor in the wrong copy or the wrong gene (LT159837).
                s, sm, sc = (p + HMM_START_SHIFT if p else 0), ("hmm" if p else "record"), ""
                e, em, ec = q, ("hmm" if q < len(x) else "record"), ""
            if e - s < 50:
                tally["too short after cut"] += 1
                continue
            flank, nc2_edits = "", ""
            if em in ("primer", "loose"):
                n2 = find_anchor_cost(x[e:e + 90], nc2r, 6)
                if n2:
                    flank, nc2_edits = x[e:e + n2[0]], n2[2]
            nc1_edits = sc if sm in ("primer", "loose") else ""
            fa.write(f">{acc}\n{x[s:e]}\n")
            rows.append((acc, sm, sc, em, ec, e - s, flank, nc1_edits, nc2_edits))
            tally[f"{sm}/{em}"] += 1
    with open(prefix + ".tsv", "w") as fh:
        fh.write("accession\tstart\tstart_edits\tend\tend_edits\tlength\tflank_28S"
                 "\tnc1_edits\tnc2_edits\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    for k, n in tally.most_common():
        print(f"{n:6d}  {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
