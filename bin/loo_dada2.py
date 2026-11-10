#!/usr/bin/env python3
"""Leave-one-out for dada2::assignTaxonomy or DECIPHER IdTaxa, on bin/its2_leave_one_out.py's queries.

    bin/loo_dada2.py <build dir> <leave-one-out prefix> <out prefix> [--folds 10] [--threads 6]
        [--method assignTaxonomy|IdTaxa]

The same queries, modes and views as <leave-one-out prefix>.tsv, so the
two classifiers' error rates compare directly. assignTaxonomy retrains on
every call, so leaving one species out at a time (1,200 trainings) is too
slow; folds instead:

  hidden   species are split into --folds groups; for each, every record
           holding a member of one of its species is removed and its
           species' queries are classified.
  present  queries are split into --folds groups, no two of one species
           in a group; for each, the group's own records are removed (the
           rest of each species stays) and the group classified.

Queries as assignTaxonomy sees reads: whole = the record; 2x150 = window 1,
10 N, window 2 (as mergePairs(justConcatenate = TRUE) writes them).
assignTaxonomy at minBoot = 80; IdTaxa as the nemabiome.ca workflow runs
it (bin/idtaxa.R: threshold 60, 100 bootstraps). Writes <out>.tsv (query, mode,
view, rank, name) and <out>.txt (rates beside the gap rule's).
"""

from __future__ import annotations

import argparse
import collections
import csv
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import RANKS, read_fasta  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def write_ref(path, keep, seqs, tax):
    with open(path, "w") as fh:
        for i in keep:
            ranks = []
            for r in RANKS:
                if not tax[i][r]:
                    break
                ranks.append(tax[i][r].replace(" ", "_"))
            fh.write(">" + ";".join(ranks) + ";\n" + seqs[i] + "\n")


def run(ref, queries, threads, d, method="assignTaxonomy"):
    q = os.path.join(d, "q.tsv")
    with open(q, "w") as fh:
        fh.write("asv\tsequence\n")
        for k, s in queries:
            fh.write(f"{k}\t{s}\n")
    out = os.path.join(d, "t.tsv")
    script = "assign_taxonomy.R" if method == "assignTaxonomy" else "idtaxa.R"
    subprocess.run(["Rscript", os.path.join(HERE, "bin", script), q, ref, out],
                   env=dict(os.environ, THREADS=str(threads)), check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    res = {}
    for r in csv.DictReader(open(out), delimiter="\t"):
        cols = [c for c in r if not c.startswith("boot") and c != "asv"][:8]
        rank, name = "", ""
        for i, c in enumerate(cols):
            if r[c] and r[c] != "NA":
                rank, name = RANKS[i], r[c].replace("_", " ")
            else:
                break
        res[r["asv"]] = (rank, name)
    return res


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("loo")
    ap.add_argument("out")
    ap.add_argument("--folds", type=int, default=10)
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--method", choices=["assignTaxonomy", "IdTaxa"], default="assignTaxonomy")
    a = ap.parse_args()
    seqs = read_fasta(os.path.join(a.build, "refdb.fasta"))
    tax = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.build, "taxonomy.tsv")),
                                               delimiter="\t")}
    holds = collections.defaultdict(set)
    for r in csv.DictReader(open(os.path.join(a.build, "members.tsv")), delimiter="\t"):
        if r["role"] == "kept" and r["id"]:
            holds[r["label"].split(";")[-1]].add(r["id"])
    loo = list(csv.DictReader(open(a.loo + ".tsv"), delimiter="\t"))
    views = sorted({r["view"] for r in loo})

    def query_seq(qi, view):
        s = seqs[qi]
        if view == "whole":
            return s
        R = int(view.split("x")[1])
        w1, w2 = R - 20, R - 67
        return s[:w1] + "N" * 10 + s[-w2:]

    rng = random.Random(11)
    species = sorted({r["species"] for r in loo if r["mode"] == "hidden"})
    rng.shuffle(species)
    sp_fold = {sp: k % a.folds for k, sp in enumerate(species)}
    q_fold = {}
    by_sp = collections.defaultdict(list)
    for r in loo:
        if r["mode"] == "present":
            by_sp[r["species"]].append(r["query"])
    for sp, qs in by_sp.items():
        start = rng.randrange(a.folds)
        for k, q in enumerate(sorted(set(qs))):
            q_fold[q] = (start + k) % a.folds

    results = {}
    with tempfile.TemporaryDirectory() as d:
        ref = os.path.join(d, "ref.fasta")
        for mode in ("hidden", "present"):
            for k in range(a.folds):
                rows = [r for r in loo if r["mode"] == mode and
                        (sp_fold.get(r["species"]) if mode == "hidden" else q_fold.get(r["query"])) == k]
                if not rows:
                    continue
                if mode == "hidden":
                    drop = set().union(*(holds[sp] | set() for sp in {r["species"] for r in rows}))
                    drop |= {r["query"] for r in rows}
                else:
                    drop = {r["query"] for r in rows}
                write_ref(ref, [i for i in seqs if i not in drop], seqs, tax)
                qs = [(f"{r['query']}|{mode}|{r['view']}", query_seq(r["query"], r["view"])) for r in rows]
                res = run(ref, qs, a.threads, d, a.method)
                for key, v in res.items():
                    results[key] = v
                print(f"{mode} fold {k + 1}/{a.folds}: {len(rows)} queries", file=sys.stderr)

    with open(a.out + ".tsv", "w") as fh:
        fh.write("query\tmode\tview\trank\tname\n")
        for key, (rank, name) in sorted(results.items()):
            q, mode, view = key.split("|")
            fh.write(f"{q}\t{mode}\t{view}\t{rank or 'unassigned'}\t{name}\n")

    def outcome(rank, name, t):
        if rank == "species":
            return "right species" if name == t["species"] else "WRONG species"
        if rank == "genus":
            return "genus, right" if name == t["genus"] else "WRONG genus"
        return "family or higher" if rank else "unassigned"
    cats = ["right species", "WRONG species", "genus, right", "WRONG genus", "family or higher",
            "unassigned"]
    lines = [f"Leave-one-out, {a.method} beside the gap rule, same queries "
             f"({a.folds} folds for assignTaxonomy).", ""]
    for mode in ("hidden", "present"):
        for view in views:
            sel = [r for r in loo if r["mode"] == mode and r["view"] == view]
            if not sel:
                continue
            lines.append(f"== {mode}, {view} ({len(sel):,} queries)")
            lines.append("   classifier     " + "".join(f"{c:>17s}" for c in cats))
            for lab, get in ((a.method, lambda r: results.get(f"{r['query']}|{mode}|{view}", ("", ""))),
                             ("gap rule", lambda r: (r["gap_rank"] if r["gap_rank"] != "unassigned" else "",
                                                     r["gap_name"]))):
                c = collections.Counter(outcome(*get(r), tax[r["query"]]) for r in sel)
                lines.append(f"   {lab:14s} " + "".join(f"{100 * c[x] / len(sel):16.1f}%" for x in cats))
            lines.append("")
    text = "\n".join(lines)
    open(a.out + ".txt", "w").write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
