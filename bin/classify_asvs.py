#!/usr/bin/env python3
"""Span-scored, tie-aware calls for DADA2 ASVs, against a build or a dada2 reference.

    bin/classify_asvs.py <asvs.tsv> <build dir | dada2 ref.fasta> <out.tsv> [--workers N]

The other arm of the comparison with assignTaxonomy (bin/dada2_nemabiome.R):
the same ASVs, scored as bin/its2_classify.py scores reads.

  - a merged ASV (primers trimmed: ITS2 plus the 28S stretch) is scored
    whole, each reference placed inside it with the ASV's extra bases free;
  - a concatenated ASV (forward + 10 N + reverse, from 2x150 pairs that do
    not overlap) is split into its two windows, scored against the start
    and the end of each reference;
  - an ASV whose best minimap2 hit is on the reverse strand (a library
    sequenced in mixed orientation) is turned round first.

A dada2 reference (headers of ;-separated ranks) is read as a build:
identical sequences are collapsed, and a record keeps only the ranks its
labels agree on - so either database can be used with either classifier.

Writes one row per ASV: edits, span, identity, ties, gap, nearest (the
tied species and their identity) and next (the nearest other species and
how many edits further it is), the call under each rule of
its2_classify.py (v4: amphid's V4 limits; exact: species
only at 0 edits, genus at <= 1; gap: >= 99% and 4 edits
clearer of any other species - the gap column is that margin; tiered:
the gap rule's species calls plus tentative ones at >= 98%, flagged in
the tentative column), and the lineage within its species that the nearest records carry
(curation/lineages.tsv), if any.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import os
import sys
import tempfile
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import RANKS, merge, read_fasta  # noqa: E402
from its2_classify import (RULES, _init, _score, call, describe, lineage,  # noqa: E402
                           load_taxonomy, minimap_candidates)
from seqmatch import revcomp  # noqa: E402

SEP = "N" * 10


def load_reference(path, tmp):
    """(refdb.fasta path, {id: sequence}, {id: taxonomy row})."""
    if os.path.isdir(path):
        tax = load_taxonomy(path)
        return os.path.join(path, "refdb.fasta"), read_fasta(os.path.join(path, "refdb.fasta")), tax
    labels = collections.defaultdict(list)
    head, buf = None, []

    def flush():
        if head is not None:
            labels["".join(buf).upper()].append(head[:len(RANKS)])
    # nemabiome's dada2 export wraps sequences at 80 bases a line.
    for line in open(path):
        if line.startswith(">"):
            flush()
            head = [x.replace("_", " ") for x in line[1:].strip().rstrip(";").split(";")]
            head += [""] * (len(RANKS) - len(head))
            buf = []
        else:
            buf.append(line.strip())
    flush()
    seqs, tax = {}, {}
    for s, labs in labels.items():
        rid = "R_" + hashlib.sha1(s.encode()).hexdigest()[:12]
        top = merge(labs)
        seqs[rid] = s
        tax[rid] = dict(zip(RANKS, top), id=rid)
    fa = os.path.join(tmp, "ref.fasta")
    with open(fa, "w") as fh:
        for rid, s in seqs.items():
            fh.write(f">{rid}\n{s}\n")
    return fa, seqs, tax


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("asvs")
    ap.add_argument("reference")
    ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    asvs = [(r["asv"], r["sequence"]) for r in csv.DictReader(open(a.asvs), delimiter="\t")]
    with tempfile.TemporaryDirectory() as tmp:
        fa, refs, tax = load_reference(a.reference, tmp)
        minus = set()
        cand = minimap_candidates(fa, [s for _, s in asvs], a.workers, minus)
    print(f"{len(minus)} of {len(asvs)} ASVs on the reverse strand, turned", file=sys.stderr)
    jobs = []
    for i, (name, s) in enumerate(asvs):
        c = cand.get(i, [])
        if i in minus:
            s = revcomp(s)
        if SEP in s:
            x, _, y = s.partition(SEP)
            jobs.append((i, "window", x, y.strip("N"), c))
        else:
            jobs.append((i, "whole", s, "", c))
    with Pool(a.workers, initializer=_init, initargs=(refs, tax)) as pool:
        res = {k: (d, sp, t, g, cv, nx) for k, d, sp, t, g, cv, nx in pool.imap_unordered(_score, jobs, chunksize=50)}
    with open(a.out, "w") as fh:
        fh.write("asv\tmode\tedits\tspan\tidentity\tties\tgap\tnearest\tnext\t"
                 + "\t".join(f"{r}_rank\t{r}_name" for r in RULES) + "\ttentative\tlineage\n")
        for i, (name, s) in enumerate(asvs):
            d, sp, t, g, cv, nx = res[i]
            near, nxt = describe(tax, t, d, sp, nx)
            ident = "" if d is None or not sp else f"{1 - d / sp:.4f}"
            calls = [call(tax, t, d, sp, r, g, cv) for r in RULES]
            by = dict(zip(RULES, calls))
            # a species named only by the tiered rule's lower floor
            tent = "yes" if by["tiered"][0] == "species" and by["gap"][0] != "species" else ""
            fh.write(f"{name}\t{'window' if SEP in s else 'whole'}\t{'' if d is None else d}\t"
                     f"{sp}\t{ident}\t{len(t)}\t{'' if g is None else g}\t{near}\t{nxt}\t"
                     + "\t".join(f"{rk or 'unassigned'}\t{nm}" for rk, nm in calls)
                     + f"\t{tent}\t{lineage(tax, t)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
