#!/usr/bin/env python3
"""Write a build as a dada2 assignTaxonomy reference.

    bin/export_dada2.py <build dir> <out.fasta>

The format of nemabiome's own dada2 export: header = the eight ranks
(superkingdom ... species), ;-separated, spaces in species names as
underscores; sequence below. A record that stops above species (its
members disagree, or it is tied) is written only to the rank it reaches,
so assignTaxonomy cannot name below it.
"""

from __future__ import annotations

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import RANKS, read_fasta  # noqa: E402


def main() -> int:
    build, out = sys.argv[1:3]
    seqs = read_fasta(os.path.join(build, "refdb.fasta"))
    n = 0
    with open(out, "w") as fh, open(os.path.join(build, "taxonomy.tsv")) as tf:
        for t in csv.DictReader(tf, delimiter="\t"):
            ranks = []
            for r in RANKS:
                if not t[r]:
                    break
                ranks.append(t[r].replace(" ", "_"))
            fh.write(">" + ";".join(ranks) + ";\n" + seqs[t["id"]] + "\n")
            n += 1
    print(f"{n} records", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
