#!/usr/bin/env python3
"""Recall, precision and over-classification for the leave-one-out arms.

    bin/loo_metrics.py <build dir> <loo.tsv> [<label>=<k-mer loo.tsv> ...]

<loo.tsv> is bin/its2_leave_one_out.py's output (the v4, exact, gap and
tiered rules); each k-mer file is bin/loo_dada2.py's (query, mode, view, rank,
name) for the same queries. Per view (whole, 2xR), at species and genus:

  present   the query's species has other records in the database
              recall     right calls / queries
              precision  right calls / calls made at that rank or deeper
  hidden    every record of the query's species removed
              species over-classification: any species call / queries
              (all are wrong)
            and at genus, split by whether the genus still has another
            species in the database:
              genus kept    genus recall and precision, as above
              genus gone    genus over-classification: any genus or
                            species call / queries

A species call is right at genus when its genus is the query's.
"""

from __future__ import annotations

import collections
import csv
import os
import sys

RULES = ("v4", "exact", "gap", "tiered")


def main() -> int:
    build, ours = sys.argv[1:3]
    kmer = [x.split("=", 1) for x in sys.argv[3:]]
    tax = list(csv.DictReader(open(os.path.join(build, "taxonomy.tsv")), delimiter="\t"))
    genus_of = {t["species"]: t["genus"] for t in tax if t["species"]}
    named = collections.defaultdict(set)
    for t in tax:
        if t["species"] and t["species_status"] == "named":
            named[t["genus"]].add(t["species"])

    def call_genus(rank, name):
        if rank == "species":
            return genus_of.get(name, name.split()[0])
        return name if rank == "genus" else ""

    truth, calls = {}, collections.defaultdict(dict)
    for r in csv.DictReader(open(ours), delimiter="\t"):
        k = (r["query"], r["mode"], r["view"])
        truth[k] = (r["species"], r["genus"])
        for rule in RULES:
            if f"{rule}_rank" in r:      # tiered: since v0.1.3
                calls[rule][k] = (r[f"{rule}_rank"], r[f"{rule}_name"])
    for label, path in kmer:
        for r in csv.DictReader(open(path), delimiter="\t"):
            k = (r["query"], r["mode"], r["view"])
            calls[label][k] = (r["rank"], r["name"])

    views = sorted({k[2] for k in truth}, key=lambda v: (v != "whole", v))
    print(f"Leave-one-out, {len({k[0] for k in truth}):,} queries, {build}.")
    print("Recall = right / queries; precision = right / calls made at that rank or deeper;")
    print("over-classification = calls made where the right answer is not in the database.\n")
    for view in views:
        keys = {m: [k for k in truth if k[1] == m and k[2] == view] for m in ("present", "hidden")}
        kept = [k for k in keys["hidden"] if named[truth[k][1]] - {truth[k][0]}]
        gone = [k for k in keys["hidden"] if k not in set(kept)]
        print(f"== {view}: present {len(keys['present']):,} queries; hidden {len(keys['hidden']):,} "
              f"(genus kept {len(kept):,}, genus gone {len(gone):,})")
        head = ["species recall", "precision", "genus recall", "precision",
                "over-class. sp", "genus kept: rec", "prec", "genus gone: over"]
        print(f"   {'arm':16s}" + "".join(f"{h:>17s}" for h in head))
        for arm, c in calls.items():
            def rp(ks, rank):
                n = right = made = 0
                for k in ks:
                    if k not in c:
                        continue
                    n += 1
                    sp, ge = truth[k]
                    rk, nm = c[k]
                    if rank == "species":
                        if rk == "species":
                            made += 1
                            right += nm == sp
                    elif rk in ("species", "genus"):
                        made += 1
                        right += call_genus(rk, nm) == ge
                return (right / n if n else 0), (right / made if made else float("nan")), made / n if n else 0
            sr, sprec, _ = rp(keys["present"], "species")
            gr, gprec, _ = rp(keys["present"], "genus")
            _, _, over_sp = rp(keys["hidden"], "species")
            kr, kprec, _ = rp(kept, "genus")
            _, _, over_g = rp(gone, "genus")
            vals = [sr, sprec, gr, gprec, over_sp, kr, kprec, over_g]
            print(f"   {arm:16s}" + "".join(f"{100 * v:16.1f}%" for v in vals))
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
