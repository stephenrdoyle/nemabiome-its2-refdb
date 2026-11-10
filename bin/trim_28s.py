#!/usr/bin/env python3
"""Cut the 28S stretch off DADA2 ASVs, for k-mer classifiers.

    bin/trim_28s.py <asvs.tsv> <out.tsv>

A primer-trimmed NC1/NC2 ASV ends in the ~47-base 28S stretch before
NC2; this database's records stop at its start (README: What a record
is). bin/classify_asvs.py scores over the reference's span, so the
stretch costs nothing there. assignTaxonomy and IdTaxa compare k-mers
over the whole ASV, so against this database every ASV carries ~47
bases no reference has, and calls are lost: run them on ASVs cut here.

The cut is at bin/its2.py's 28S anchor (<= 12 edits, within the last 90
bases), taken only if it removes 35-60 bases; else the last 47 bases go.
A concatenated ASV (forward + 10 N + reverse, from 2x150 pairs) is cut
in its reverse part. ASVs with under 100 bases left are not ITS2 and are
dropped (a few short junk ASVs; they also overflow IdTaxa's memory).
Writes asv, sequence.
"""

from __future__ import annotations

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from its2 import S28  # noqa: E402
from seqmatch import find_anchor_cost  # noqa: E402

STRETCH, LOOK, MAX_EDITS = 47, 90, 12
CUT_RANGE = (35, 60)
MIN_LEN = 100
SEP = "N" * 10


def trim(s):
    """(trimmed sequence, whether the anchor was found)."""
    head, sep, tail = s.partition(SEP)
    part = tail if sep else head
    hit = find_anchor_cost(part[-LOOK:], S28, MAX_EDITS) if len(part) >= LOOK else None
    cut = len(part) - LOOK + hit[0] if hit else None
    found = cut is not None and CUT_RANGE[0] <= len(part) - cut <= CUT_RANGE[1]
    if not found:
        cut = len(part) - STRETCH       # anchor missing or out of place
    part = part[:cut]
    return (head + sep + part if sep else part), found


def main() -> int:
    src, out = sys.argv[1:3]
    n = found = dropped = 0
    with open(out, "w") as fh:
        fh.write("asv\tsequence\n")
        for r in csv.DictReader(open(src), delimiter="\t"):
            n += 1
            s, f = trim(r["sequence"])
            found += f
            if len(s.replace("N", "")) < MIN_LEN:
                dropped += 1
                continue
            fh.write(f"{r['asv']}\t{s}\n")
    print(f"{src}: {n} ASVs, 28S anchor found in {found}, {dropped} under {MIN_LEN} bases dropped",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
