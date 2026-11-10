#!/usr/bin/env python3
"""Where the k-mer classifiers' extra species calls come from.

    bin/loo_kmer_vs_gap.py <loo.tsv> <assignTaxonomy loo.tsv> <IdTaxa loo.tsv> > <out.txt>

From bin/its2_leave_one_out.py's and bin/loo_dada2.py's outputs for the
same queries (whole reads). Two tables:

  1. Present queries, by why the gap rule stops short of species, with
     how often each k-mer classifier names the right or wrong species on
     those reads (points of all present queries).
  2. By identity of the best hit, the k-mer classifiers' species-call
     rate on present queries whose own species is best (right calls) and
     on hidden queries (wrong calls): does a k-mer classifier tell a
     divergent variant from a missing species better than identity does?
"""
import collections, csv, sys

ours = {(r["query"], r["mode"], r["view"]): r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
km = {}
for lab, p in (("assignTaxonomy", sys.argv[2]), ("IdTaxa", sys.argv[3])):
    km[lab] = {(r["query"], r["mode"], r["view"]): (r["rank"], r["name"])
               for r in csv.DictReader(open(p), delimiter="\t")}


def reason(r):
    if r["gap_rank"] == "species":
        return "gap rule names species"
    if not r["identity"]:
        return "no hit"
    ties = r["tied_species"].split(";")
    if len(ties) > 1:
        return "tied with another species"
    if r["species"] not in ties:
        return "another species closest"
    if float(r["identity"]) < 0.99:
        return "own species best, < 99%"
    return "own >= 99%, sister 0-3 edits off"


ORDER = ["gap rule names species", "own species best, < 99%", "own >= 99%, sister 0-3 edits off",
         "another species closest", "tied with another species", "no hit"]
print("Whole reads. Points of all queries of the mode; right / wrong species calls.\n")
for mode in ("present", "hidden"):
    ks = [k for k in ours if k[1] == mode and k[2] == "whole"]
    t = collections.defaultdict(collections.Counter)
    for k in ks:
        r = ours[k]
        rs = reason(r)
        t[rs]["n"] += 1
        for lab in km:
            rk, nm = km[lab].get(k, ("", ""))
            if rk == "species":
                t[rs][lab, nm == r["species"]] += 1
    n = len(ks)
    print(f"== {mode} ({n:,} queries)")
    print(f"   {'gap rule outcome':36s} {'share':>6s}" + "".join(f"  {lab:>20s}" for lab in km))
    for rs in ORDER:
        if rs in t:
            c = t[rs]
            print(f"   {rs:36s} {100 * c['n'] / n:5.1f}%" + "".join(
                f"  {100 * c[lab, True] / n:9.1f} / {100 * c[lab, False] / n:5.1f}" for lab in km))
    print()

print("Species-call rate by identity of the best hit: present (own species best; calls right) /")
print("hidden (species missing; calls wrong).")
print(f"   {'identity':9s} {'present':>8s} {'hidden':>7s}" + "".join(f"  {lab:>20s}" for lab in km))
for lo, hi, lab in ((0.99, 1.01, ">= 99%"), (0.98, 0.99, "98-99%"), (0.97, 0.98, "97-98%"),
                    (0.95, 0.97, "95-97%"), (0.90, 0.95, "90-95%")):
    P = [k for k, r in ours.items() if k[2] == "whole" and r["mode"] == "present" and r["identity"]
         and lo <= float(r["identity"]) < hi and r["species"] in r["tied_species"].split(";")]
    H = [k for k, r in ours.items() if k[2] == "whole" and r["mode"] == "hidden" and r["identity"]
         and lo <= float(r["identity"]) < hi]
    row = f"   {lab:9s} {len(P):8,d} {len(H):7,d}"
    for m in km:
        p = sum(km[m].get(k, ("",))[0] == "species" for k in P)
        h = sum(km[m].get(k, ("",))[0] == "species" for k in H)
        row += f"  {100 * p / len(P):6.0f}% / {100 * h / len(H):3.0f}% ({p}:{h})"
    print(row)
