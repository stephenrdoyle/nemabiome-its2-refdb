#!/usr/bin/env python3
"""Leave-one-out: what the classifier calls a read whose species is missing.

    bin/its2_leave_one_out.py <build dir> <out prefix> [--per-species 5] [--read-lengths 150] [--workers N]

As vespa-18s-refdb's bin/simulate_leave_one_out.py, for the ITS2 classifier
(bin/its2_classify.py). Queries are full-length (complete or
near_complete) species-level records, up to --per-species per species,
each read as a perfect read: whole (nanopore, or merged pairs), and
through the 2xR windows (first R-20 bases; last R-67 plus the 28S
stretch, here N).

  hidden   every record holding a member of the query's species is removed
           (shared records included): any species named is wrong; the
           right answer is the genus or above.
  present  only the query's own record is removed, for species with at
           least one other record: the species is still there.

Writes <prefix>.tsv: per query, mode, view; the best identity, gap, tied
species and their common ancestor - the raw material, so call rules can
be compared without rescoring (bin/its2_classify.py's rules are applied
in <prefix>.txt).
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import random
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import RANKS, deepest, merge_records, read_fasta  # noqa: E402
from its2_classify import RULES, STRETCH_28S, _init, _score, call, minimap_candidates  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("prefix")
    ap.add_argument("--per-species", type=int, default=5)
    ap.add_argument("--read-lengths", default="150")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    seqs = read_fasta(os.path.join(a.build, "refdb.fasta"))
    tax = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.build, "taxonomy.tsv")),
                                               delimiter="\t")}
    holds = collections.defaultdict(set)      # species -> records holding a member
    for r in csv.DictReader(open(os.path.join(a.build, "members.tsv")), delimiter="\t"):
        if r["role"] == "kept" and r["id"]:
            holds[r["label"].split(";")[-1]].add(r["id"])
    by_sp = collections.defaultdict(list)
    for i, t in tax.items():
        if t["deepest_rank"] == "species" and t["grade"] in ("complete", "near_complete"):
            by_sp[t["species"]].append(i)
    rng = random.Random(7)
    queries = [q for sp in sorted(by_sp) for q in
               rng.sample(sorted(by_sp[sp]), min(a.per_species, len(by_sp[sp])))]
    print(f"{len(queries):,} queries from {len(by_sp):,} species", file=sys.stderr)

    found = minimap_candidates(os.path.join(a.build, "refdb.fasta"), [seqs[i] for i in queries], a.workers)
    cand = {qi: found.get(n, []) for n, qi in enumerate(queries)}

    views = [("whole", None)] + [(f"2x{r}", int(r)) for r in a.read_lengths.split(",")]
    jobs, meta = [], []
    for qi in queries:
        sp = tax[qi]["species"]
        for mode, drop in (("hidden", holds[sp] | {qi}), ("present", {qi})):
            if mode == "present" and not (holds[sp] - {qi}):
                continue      # no other record of the species: present = hidden
            cs = [c for c in cand[qi] if c not in drop]
            s = seqs[qi]
            for view, r in views:
                if r is None:
                    jobs.append(((qi, mode, view), "whole", s + "N" * STRETCH_28S, "", cs))
                else:
                    w1, w2 = r - 20, r - 20 - STRETCH_28S
                    if len(s) <= w1 + w2:
                        continue      # the windows see it whole
                    jobs.append(((qi, mode, view), "window", s[:w1], s[-w2:] + "N" * STRETCH_28S, cs))
    with Pool(a.workers, initializer=_init, initargs=(seqs, tax)) as pool:
        res = list(pool.imap_unordered(_score, jobs, chunksize=50))

    rows = []
    with open(a.prefix + ".tsv", "w") as fh:
        fh.write("query\tspecies\tgenus\tfamily\tmode\tview\tedits\tspan\tidentity\tgap\t"
                 "tied_species\tlca_rank\tlca_name\t" + "\t".join(f"{r}_rank\t{r}_name" for r in RULES)
                 + "\n")
        for (qi, mode, view), dd, sp, ties, gap, cv, _ in res:
            t = tax[qi]
            lab = merge_records([[tax[x][k] for k in RANKS] for x in ties]) if ties else [""] * 8
            lr, ln = deepest(lab)
            calls = [call(tax, ties, dd, sp, r, gap, cv) for r in RULES]
            ident = "" if dd is None or not sp else f"{1 - dd / sp:.4f}"
            row = [qi, t["species"], t["genus"], t["family"], mode, view, "" if dd is None else dd,
                   sp, ident, "" if gap is None else gap,
                   ";".join(sorted({tax[x]["species"] for x in ties})), lr, ln]
            row += [x for rk, nm in calls for x in (rk or "unassigned", nm)]
            rows.append(dict(zip(["query", "species", "genus", "family", "mode", "view"], row[:6]),
                             calls={r: c for r, c in zip(RULES, calls)}))
            fh.write("\t".join(map(str, row)) + "\n")

    def outcome(r, rule):
        rank, name = r["calls"][rule]
        t = tax[r["query"]]
        if rank == "species":
            return "right species" if name == t["species"] else "WRONG species"
        if rank == "genus":
            return "genus, right" if name == t["genus"] else "WRONG genus"
        return "family or higher" if rank else "unassigned"
    cats = ["right species", "WRONG species", "genus, right", "WRONG genus", "family or higher",
            "unassigned"]
    lines = [f"Leave-one-out, {len(queries):,} queries ({a.per_species} per species at most) "
             f"from {len(by_sp):,} species in {a.build}.", ""]
    for mode in ("hidden", "present"):
        for view, _ in views:
            sel = [r for r in rows if r["mode"] == mode and r["view"] == view]
            if not sel:
                continue
            lines.append(f"== {mode}, {view} ({len(sel):,} queries)")
            lines.append("   rule   " + "".join(f"{c:>17s}" for c in cats))
            for rule in RULES:
                c = collections.Counter(outcome(r, rule) for r in sel)
                lines.append(f"   {rule:6s} " + "".join(f"{100 * c[k] / len(sel):16.1f}%" for k in cats))
            lines.append("")
    text = "\n".join(lines)
    open(a.prefix + ".txt", "w").write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
