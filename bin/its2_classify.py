#!/usr/bin/env python3
"""Classify NC1/NC2 read pairs against a build, to test the database.

    bin/its2_classify.py <build dir> <out prefix> <R1.fastq.gz> [<R1> ...]
        [--read-lengths 150,250] [--workers N]

Not amphid: a small stand-in that scores reads the way
docs/feedback/its2_notes_for_amphid.md asks amphid to, so the database
can be judged on real data before the pipeline changes.

Per read pair (R2 found by replacing _1. with _2.):
  - NC1 at the start of R1, NC2 at the start of R2, each after 0-3
    spacer bases and within 2 edits; pairs without both are skipped.
  - whole: R1 and R2 inserts are merged (best overlap >= 20 bases, <= 5%
    mismatches) into the amplicon - ITS2 plus the ~47-base 28S stretch.
    Each reference is placed whole inside it, the read's bases beyond
    it free (score over the reference's span).
  - 2xR, for each read length R in --read-lengths (as in vespa-18s-refdb,
    the read length is set, not assumed; default 150): the first R bases
    of each read, primers trimmed. Window 1 is R1's insert, window 2 is
    R2's insert reverse-complemented. Window 1 is scored against the
    start of each reference, window 2 against its end, the 28S stretch
    in window 2 free. Reads longer than R are cut to R, so 2x250 data
    stands in for 2x150.

Distances are exact (bit-parallel edit distance) among candidates
minimap2 proposes. Every reference at the best distance is a tie; the
call is their common ancestor, then cut back by the identity rules:
  v4      amphid's V4 limits: species >= 99.7%, genus >= 99.5%,
          family >= 95% (Illumina)
  exact   species only at 0 edits, genus <= 1, family >= 95%

Writes <prefix>.calls.tsv (run, mode, rule, rank, name, reads) and
<prefix>.reads.tsv (per unique sequence: reads, best edits, span,
identity, ties, tied names). Unique sequences are scored once.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import RANKS, deepest, infix_distance, merge_records, read_fasta  # noqa: E402
from its2 import NC1, NC2  # noqa: E402
from seqmatch import find_anchor_cost, revcomp  # noqa: E402

SPACER = 3
STRETCH_28S = 47
TOP = 40
RULES = {"v4": [("species", 0.997), ("genus", 0.995), ("family", 0.95)],
         "exact": [("species", 1.0), ("genus", None), ("family", 0.95)]}


def fastq_pairs(r1):
    r2 = r1.replace("_1.fastq", "_2.fastq")
    with gzip.open(r1, "rt") as a, gzip.open(r2, "rt") as b:
        while True:
            h = a.readline()
            if not h:
                return
            s1 = a.readline().strip()
            a.readline()
            a.readline()
            b.readline()
            s2 = b.readline().strip()
            b.readline()
            b.readline()
            yield s1, s2


def trim(seq, primer):
    h = find_anchor_cost(seq[:len(primer) + SPACER + 2], primer, 2)
    return seq[h[1]:] if h else None


def merge(i1, i2rc):
    """Merge R1 insert with reverse-complemented R2 insert."""
    best = None
    for ov in range(min(len(i1), len(i2rc)), 19, -1):
        a, b = i1[-ov:], i2rc[:ov]
        mm = sum(x != y for x, y in zip(a, b))
        if mm <= 0.05 * ov:
            best = ov
            break
    return i1 + i2rc[best:] if best else None


_R = _T = None


def _init(refs, tax):
    global _R, _T
    _R, _T = refs, tax


def _score(job):
    """(key, best edits, span, ties) for one unique read or window pair."""
    key, mode, a, b, cands = job
    scores = []
    for c in cands:
        r = _R[c]
        if mode == "whole":
            d, span = infix_distance(r, a), len(r)
        else:
            k1 = min(len(a), len(r))
            k2 = max(0, min(len(b) - STRETCH_28S, len(r)))
            d1 = min(infix_distance(r[:k1], a), infix_distance(a[:k1], r[:k1 + 10]))
            d2 = infix_distance(r[len(r) - k2:], b) if k2 else 0
            d, span = d1 + d2, k1 + k2
        scores.append((d, span, c))
    if not scores:
        return key, None, 0, []
    d0 = min(s[0] for s in scores)
    ties = [s for s in scores if s[0] == d0]
    return key, d0, max(s[1] for s in ties), [s[2] for s in ties]


def call(tax, ties, d, span, rule):
    lab = merge_records([[tax[t][k] for k in RANKS] for t in ties])
    if d is None or span == 0:
        return "", "unassigned"
    ident = 1 - d / span
    allowed = None
    for rank, lim in RULES[rule]:
        ok = (d <= 1) if lim is None else ident >= lim
        if ok:
            allowed = rank
            break
    if not allowed:
        return "", "unassigned"
    cut = RANKS.index(allowed)
    lab = [x if i <= cut else "" for i, x in enumerate(lab)]
    return deepest(lab)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("prefix")
    ap.add_argument("r1", nargs="+")
    ap.add_argument("--read-lengths", default="150",
                    help="comma-separated Illumina read lengths to test (default 150)")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    lengths = [int(x) for x in a.read_lengths.split(",")]
    modes = [f"2x{r}" for r in lengths]
    refs = read_fasta(os.path.join(a.build, "refdb.fasta"))
    tax = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.build, "taxonomy.tsv")),
                                               delimiter="\t")}

    # --- reads -> unique whole amplicons and unique 2x150 window pairs
    per_run = {}
    uniq = {m: collections.Counter() for m in ["whole"] + modes}
    stats = collections.Counter()
    for r1 in a.r1:
        run = os.path.basename(r1).split("_")[0]
        cnt = {m: collections.Counter() for m in uniq}
        for s1, s2 in fastq_pairs(r1):
            stats["pairs"] += 1
            i1, i2 = trim(s1, NC1), trim(s2, NC2)
            if not i1 or not i2 or "N" in s1 or "N" in s2:
                continue
            stats["primed"] += 1
            for r, m in zip(lengths, modes):
                cnt[m][(i1[:r - (len(s1) - len(i1))], revcomp(i2[:r - (len(s2) - len(i2))]))] += 1
            j = merge(i1, revcomp(i2))
            if j:
                stats["merged"] += 1
                cnt["whole"][j] += 1
        per_run[run] = cnt
        for m in uniq:
            uniq[m].update(cnt[m])
    print(f"{stats['pairs']:,} pairs; {stats['primed']:,} with both primers; "
          f"{stats['merged']:,} merged; unique: "
          + ", ".join(f"{m} {len(c):,}" for m, c in uniq.items()), file=sys.stderr)

    # --- candidates from minimap2 (whole amplicon, or the two windows joined)
    def candidates(items):
        with tempfile.TemporaryDirectory() as d:
            q = os.path.join(d, "q.fa")
            with open(q, "w") as fh:
                for i, s in enumerate(items):
                    fh.write(f">{i}\n{s}\n")
            out = subprocess.run(["minimap2", "-t", str(a.workers), "-k", "11", "-w", "5",
                                  "-f", "1000000", "-N", str(TOP), "-p", "0.5",
                                  os.path.join(a.build, "refdb.fasta"), q],
                                 capture_output=True, text=True, check=True).stdout
        c = collections.defaultdict(dict)
        for line in out.splitlines():
            f = line.split("\t")
            c[int(f[0])][f[5]] = max(c[int(f[0])].get(f[5], 0), int(f[9]))
        return {i: sorted(v, key=lambda t: -v[t])[:TOP] for i, v in c.items()}

    items = {m: list(c) for m, c in uniq.items()}
    jobs = []
    for m, its in items.items():
        cand = candidates(its if m == "whole" else [x + "N" * 20 + y for x, y in its])
        for i, s in enumerate(its):
            x, y = (s, "") if m == "whole" else s
            jobs.append(((m, i), "whole" if m == "whole" else "window", x, y, cand.get(i, [])))
    with Pool(a.workers, initializer=_init, initargs=(refs, tax)) as pool:
        res = {k: (d, sp, t) for k, d, sp, t in pool.imap_unordered(_score, jobs, chunksize=200)}

    idx = {m: {s: i for i, s in enumerate(its)} for m, its in items.items()}
    with open(a.prefix + ".reads.tsv", "w") as fh:
        fh.write("mode\tkey\treads\tedits\tspan\tidentity\tties\ttied_species\n")
        for mode, its in items.items():
            for i, s in enumerate(its):
                d, sp, t = res[(mode, i)]
                names = sorted({tax[x]["species"] or tax[x]["deepest_rank"] for x in t})
                ident = "" if d is None or not sp else f"{1 - d / sp:.4f}"
                fh.write(f"{mode}\t{i}\t{uniq[mode][s]}\t{'' if d is None else d}\t{sp}\t{ident}\t{len(t)}\t"
                         f"{';'.join(names[:5])}\n")
    rows = collections.Counter()
    for run, cnt in per_run.items():
        for mode, c in cnt.items():
            for s, n in c.items():
                d, sp, t = res[(mode, idx[mode][s])]
                for rule in RULES:
                    rank, name = call(tax, t, d, sp, rule)
                    rows[(run, mode, rule, rank or "unassigned", name)] += n
    with open(a.prefix + ".calls.tsv", "w") as fh:
        fh.write("run\tmode\trule\trank\tname\treads\n")
        for k, n in sorted(rows.items()):
            fh.write("\t".join(k) + f"\t{n}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
