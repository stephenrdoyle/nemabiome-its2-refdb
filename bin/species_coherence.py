#!/usr/bin/env python3
"""Records that do not sit with their own species.

    bin/species_coherence.py <build dir> [--workers N] > <out.tsv>

For every named record of a species with two or more records: `own`,
the edits to the nearest other record of its species, and `other`, the
edits to the nearest record of another species (species_margins.tsv),
both over the shorter span (bin/build.py's span_distance). A record
that sits far from its own species is either a mislabel, a partial that
ends early, or real variation the species' other records do not sample;
`own_spread` (the species' median nearest-own distance) says which is
plausible. Calls:

  displaced   other <= 3 and own >= 10: next to another species, far
              from its own - the likely mislabel, named by `other_species`
  stray       own >= 20, and >= 20 from the species' best-supported
              record (the one with most accessions, >= 3 of them, its
              `core`), and nothing within 10: far from its own species
              and from everything else - an unknown, or another species
              not in the database
  closer      other < own (and not displaced or stray): nearer another
              species than its own - check

Records of single-record species are not judged here (no own species to
compare); bin/outlier_records.py covers those that crowd another species.
Records from genome assemblies and long reads (a layer: `sources`
neither the backbone nor bin/genbank_update.py's) are rDNA copies within one worm, not labels to doubt: they
are listed as "layer" whatever their distances.
"""
from __future__ import annotations

import argparse
import collections
import csv
import os
import statistics
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import read_fasta, span_distance  # noqa: E402

_S = None


def _init(seqs):
    global _S
    _S = seqs


def _own(job):
    i, others = job
    return i, min((span_distance(_S[i], _S[j]), j) for j in others)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    b = a.build
    seqs = read_fasta(os.path.join(b, "refdb.fasta"))
    tax = {r["id"]: r for r in csv.DictReader(open(os.path.join(b, "taxonomy.tsv")), delimiter="\t")}
    mg = {r["id"]: r for r in csv.DictReader(open(os.path.join(b, "species_margins.tsv")),
                                             delimiter="\t")}
    acc = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(b, "members.tsv")), delimiter="\t"):
        if r["id"] and r["role"] != "dropped":
            acc[r["id"]].append(r["accession"])
    by_sp = collections.defaultdict(list)
    backbone = collections.Counter(t["sources"] for t in tax.values()).most_common(1)[0][0]
    # GenBank records (the backbone, and bin/genbank_update.py's) carry
    # labels to doubt; layer copies do not.
    genbank = lambda srcs: any(x == backbone or x.startswith("genbank_update") for x in srcs.split(";"))
    for i, t in tax.items():
        if t["species"] and t["species_status"] == "named":
            by_sp[t["species"]].append(i)
    jobs = [(i, [j for j in ids if j != i]) for ids in by_sp.values() if len(ids) > 1 for i in ids]
    with Pool(a.workers, initializer=_init, initargs=(seqs,)) as pool:
        own = dict(pool.imap_unordered(_own, jobs, chunksize=20))
    spread = {sp: statistics.median(own[i][0] for i in ids)
              for sp, ids in by_sp.items() if len(ids) > 1}
    core = {}
    for sp, ids in by_sp.items():
        c = max(ids, key=lambda i: (len(acc[i]), i))
        if len(acc[c]) >= 3:
            core[sp] = c
    cols = ["call", "id", "species", "n_species", "n_acc", "own", "own_spread", "to_core",
            "core_acc", "other", "other_species", "other_record_own", "grade", "accessions",
            "own_nearest_accessions"]
    print("\t".join(cols))
    rows = []
    for i, (d_own, j) in own.items():
        sp = tax[i]["species"]
        m = mg.get(i, {})
        d_oth = int(m["margin_whole"]) if m.get("margin_whole", "") != "" else None
        n = m.get("nearest_whole", "")
        o_sp = (tax[n]["species"] or f"({tax[n]['genus']} {tax[n]['species_status']})") if n else ""
        c = core.get(sp)
        d_core = "" if c is None or c == i else span_distance(seqs[i], seqs[c])
        if not genbank(tax[i]["sources"]):
            if d_own < 10 and (d_oth is None or d_oth >= d_own):
                continue
            call = "layer"
        elif d_oth is not None and d_oth <= 3 and d_own >= 10:
            call = "displaced"
        elif d_own >= 20 and d_core != "" and d_core >= 20 and (d_oth is None or d_oth > 10):
            call = "stray"
        elif d_oth is not None and d_oth < d_own:
            call = "closer"
        else:
            continue
        rows.append([call, i, sp, len(by_sp[sp]), len(acc[i]), d_own, f"{spread[sp]:g}", d_core,
                     len(acc[c]) if c else "",
                     "" if d_oth is None else d_oth, o_sp,
                     own[n][0] if n in own else "", tax[i]["grade"],
                     ",".join(acc[i][:6]), ",".join(acc[j][:3])])
    order = {"displaced": 0, "stray": 1, "closer": 2, "layer": 3}
    for r in sorted(rows, key=lambda r: (order[r[0]], r[9] if r[9] != "" else 99, -r[5])):
        print("\t".join(map(str, r)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
