#!/usr/bin/env python3
"""Accessions (no version) whose label is an informal "Genus sp." name.

    bin/informal_accessions.py <taxonomy.txt> > informal.txt

For bin/its2_partials.py, which never takes these records as a partial's
relative (bin/build.py, INFORMAL NAMES).
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import RANKS, informal  # noqa: E402

for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    if informal([r[k] for k in RANKS]):
        print(r["accn"].split(".")[0])
