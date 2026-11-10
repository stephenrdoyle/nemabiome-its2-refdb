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

Within each run the reads of each mode are denoised first (UNOISE-style,
see denoise()): a rare sequence 1-2 edits from a much commoner one is
taken as its sequencing error and folded into it, and sequences with a
single read in the run are dropped as noise. Only the centroids are scored.

Distances are exact (bit-parallel edit distance) against every
reference minimap2 proposes (minimap_candidates(): no cap). References
are ranked by identity over the span
scored (at least MIN_SPAN bases); every reference at the best identity
whose span is within SPAN_SLACK of the longest such span is a tie (a
partial matching only over a shorter span does not tie with a longer
match); the
call is their common ancestor, then cut back by the identity rules:
  v4      amphid's V4 limits: species >= 99.7%, genus >= 99.5%,
          family >= 95% (Illumina)
  exact   species only at 0 edits, genus <= 1, family >= 95%
  gap     species when >= 99% identical AND the nearest reference of any
          other species is >= 4 edits further away (exact matches too:
          sister species share variants); genus >= 97%, family >= 95%.
          The gap is reported with the call as its confidence.

Writes <prefix>.calls.tsv (run, mode, rule, rank, name, reads) and
<prefix>.reads.tsv (per unique sequence: reads, best edits, span,
identity, ties, nearest and next - see describe()). Unique sequences are scored once.
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
ALPHA = 2           # UNOISE3's default
MIN_READS = 2       # a centroid seen once in a run is noise
MIN_SPAN = 100      # bases a reference must cover to count
SPAN_SLACK = 5      # spans this close are equally strong evidence
TRIM, TRIM_STEP, SLACK = 30, 3, 10   # window overlap with partial references
RULES = {"v4": [("species", 0.997), ("genus", 0.995), ("family", 0.95)],
         "exact": [("species", 1.0), ("genus", None), ("family", 0.95)],
         "gap": [("species", "gap"), ("genus", 0.97), ("family", 0.95)]}
# Set by leave-one-out (bin/its2_leave_one_out.py, DESIGN 5.6). A species
# call needs >= 99% identity and the nearest other species >= 4 edits
# further away - exact matches included: sister species share variants,
# so an exact match with another species close by was right only 75% of
# the time. With the species in the database, 98.3% (whole) and 98.7%
# (2x150) of species calls are right; with it missing, 6.8% and 9.1% of
# its reads are named as a sister species, which no per-read rule can
# prevent. The gap is reported with each call as its confidence.
GAP_IDENTITY, GAP_EDITS = 0.99, 4


# Every reference minimap2 proposes is scored. A cap on the top N by
# matching bases hid the nearest other species from reads of species with
# hundreds of records: that neighbour is often a partial with fewer
# matches (a T. circumcincta ASV, next O. ostertagi +3, read as +9 under a
# cap of 40; issue #1). -N above the record count, so minimap2 drops
# nothing either; -p 0.3 keeps hits down to 30% of the best chain.
MM2 = ["-k", "11", "-w", "5", "-f", "1000000", "-N", "20000", "-p", "0.3"]


def minimap_candidates(fasta, queries, workers):
    """{query index: [reference ids]} for a list of query sequences."""
    with tempfile.TemporaryDirectory() as d:
        q = os.path.join(d, "q.fa")
        with open(q, "w") as fh:
            fh.writelines(f">{i}\n{s}\n" for i, s in enumerate(queries))
        out = subprocess.run(["minimap2", "-t", str(workers)] + MM2 + [fasta, q],
                             capture_output=True, text=True, check=True).stdout
    c = collections.defaultdict(set)
    for line in out.splitlines():
        f = line.split("\t", 6)
        c[int(f[0])].add(f[5])
    return {i: sorted(v) for i, v in c.items()}


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


def _key(x):
    return x if isinstance(x, str) else x[0] + "|" + x[1]


def distance(a, b, cap=3):
    """Global edit distance, near enough, capped: equal lengths are
    compared base by base first (substitutions are most Illumina errors),
    otherwise each is placed inside the other."""
    if len(a) == len(b):
        mm = 0
        for x, y in zip(a, b):
            if x != y:
                mm += 1
                if mm >= cap:
                    break
        if mm < cap:
            return mm
    return max(infix_distance(a, b), infix_distance(b, a))


def denoise(counts, alpha=ALPHA, max_d=2, min_reads=MIN_READS):
    """UNOISE-style (Edgar 2016): in order of abundance, a sequence within
    max_d edits of a more abundant centroid with at least
    2^(alpha*d + 1) times its reads is that centroid's error and is
    folded into it; otherwise it is a centroid itself. Centroids with
    fewer than min_reads reads are noise, and so are sequences with
    fewer than min_reads reads from the start: folded or not, they are
    never scored, and comparing each of them (most of the uniques) with
    every centroid is what makes denoising slow. Returns (centroid
    counts, reads folded, reads dropped as noise)."""
    cents, by_len = {}, collections.defaultdict(list)
    folded = noise = 0
    for x, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        if n < min_reads:
            noise += n
            continue
        k = _key(x)
        parent = None
        for L in range(len(k) - max_d, len(k) + max_d + 1):
            for c in by_len.get(L, ()):
                if cents[c] < n * 2 ** (alpha + 1):
                    break          # sorted by abundance: no later one is big enough
                d = distance(k, _key(c))
                if d <= max_d and cents[c] >= n * 2 ** (alpha * d + 1):
                    parent = c
                    break
            if parent is not None:
                break
        if parent is None:
            cents[x] = n
            by_len[len(k)].append(x)
        else:
            cents[parent] += n
            folded += n
    kept = {x: n for x, n in cents.items() if n >= min_reads}
    return kept, folded, noise + sum(cents.values()) - sum(kept.values())


_R = _T = None


def _init(refs, tax):
    global _R, _T
    _R, _T = refs, tax


def _window(r, w, start):
    """Edits and span of read window w against the start (or end) of
    reference r. The reference may run past the window (it is placed with
    slack), or stop short of it - a partial reference missing the ITS2's
    first or last bases: then the window's outer bases have nothing to
    meet, so up to TRIM of them are left out and only the overlap is
    scored. Without this, a partial O. radiatum reference charged 3
    edits to reads that match it exactly (a cattle set, 2x150)."""
    if not w:
        return 0, 0
    best = None
    for t in range(0, max(0, min(TRIM, len(w) - MIN_SPAN // 2)) + 1, TRIM_STEP):
        x = w[t:] if start else w[:len(w) - t]
        tail = r[:len(x) + SLACK] if start else r[-(len(x) + SLACK):]
        d = infix_distance(x, tail)
        if best is None or (1 - d / len(x)) > (1 - best[0] / best[1]) + 1e-9:
            best = (d, len(x))
        if d == 0:
            break
    return best


def _score(job):
    """(key, best edits, span, ties, gap, cover) for one unique read or
    window pair: gap is the edits by which the nearest other-species
    reference falls behind the best (None when the ties do not name one
    species); cover is the share of the read's scorable length the best
    span covers."""
    key, mode, a, b, cands = job
    scores = []
    for c in cands:
        r = _R[c]
        if mode == "whole":
            d, span = infix_distance(r, a), len(r)
        else:
            d1, k1 = _window(r, a, start=True)
            d2, k2 = _window(r, b[:max(0, len(b) - STRETCH_28S)], start=False)
            d, span = d1 + d2, k1 + k2
        scores.append((d, span, c))
    # Rank by identity over the span scored, not by raw edits: a short
    # partial reference has fewer bases to get wrong, and on edits alone a
    # 77-base fragment 32 edits off beat full-length references (every
    # 2x150 ASV against the uncut nemabiome backbone). Spans under MIN_SPAN
    # are not evidence.
    scores = [x for x in scores if x[1] >= MIN_SPAN]
    # The read's scorable length: whole, less its 28S stretch; windowed,
    # both windows less window 2's stretch.
    qlen = (len(a) - STRETCH_28S) if mode == "whole" else (len(a) + max(0, len(b) - STRETCH_28S))
    if not scores:
        return key, None, 0, [], None, 0.0, None
    best = max(1 - d / sp for d, sp, _ in scores)
    ties = [x for x in scores if 1 - x[0] / x[1] >= best - 1e-9]
    # Among equally identical references, the one covering more of the read
    # is the stronger evidence: a partial that matches only over its shorter
    # span does not tie with a full-length match (an H. placei ASV, 23% of
    # one cattle set, tied with a partial H. contortus record and stopped
    # at genus). SPAN_SLACK keeps boundary jitter from breaking real ties.
    longest = max(x[1] for x in ties)
    ties = [x for x in ties if x[1] >= longest - SPAN_SLACK]
    d0, sp0 = min((x[0], -x[1]) for x in ties)
    # The margin: how many more edits (over the same span) the nearest
    # reference NOT of the tied species needs.
    tied = {_T[x[2]]["species"] for x in ties}
    # A record naming no species counts as another species only when its
    # members name two (shared); hybrids and other unnamed records say
    # nothing about which species a read is from.
    other = [(1 - d / sp, c) for d, sp, c in scores if _T[c]["species"] not in tied
             and _T[c].get("species_status", "shared") != "unnamed"]
    nxt = None
    if other:
        oi, oc = max(other)
        nxt = (label_of(_T[oc]), round((1 - oi) * -sp0) - d0)
    gap = None
    if len(tied) == 1 and "" not in tied:
        gap = nxt[1] if nxt else -sp0 - d0
    return key, d0, -sp0, [x[2] for x in ties], gap, min(1.0, -sp0 / max(1, qlen)), nxt


def label_of(t):
    """A record's species; for a shared record the species its members
    name (Cooperia shared: oncophora / spatulata); else its deepest rank."""
    if t["species"]:
        return t["species"]
    rank, name = deepest([t[k] for k in RANKS])
    if t.get("named"):
        return f"{name} shared: " + " / ".join(s.split(" ", 1)[-1] for s in t["named"])
    return f"{name} ({rank})" if name else "unassigned"


def load_taxonomy(build):
    """taxonomy.tsv, with each shared record's member species as "named"."""
    tax = {r["id"]: r for r in csv.DictReader(open(os.path.join(build, "taxonomy.tsv")),
                                               delimiter="\t")}
    for r in csv.DictReader(open(os.path.join(build, "members.tsv")), delimiter="\t"):
        t = tax.get(r["id"])
        if t and r["role"] == "kept" and t["species_status"] == "shared":
            sp = r["label"].split(";")[-1]
            if sp and sp not in t.setdefault("named", []):
                t["named"].append(sp)
                t["named"].sort()
    return tax


def describe(tax, ties, d, span, nxt):
    """The nearest and next columns: the tied species and their identity
    (Haemonchus contortus 99.6%), and the nearest reference of any other
    species with how many edits further it is (Haemonchus placei +2) -
    so a call stopped at genus still says how close it came."""
    if d is None or not span:
        return "", ""
    # A tied shared record stands for the species its members name.
    names = sorted({n for t in ties for n in (tax[t].get("named") or [label_of(tax[t])])})
    near = " / ".join(names[:3]) + (f" / +{len(names) - 3} more" if len(names) > 3 else "")
    # No other species among the candidates minimap2 proposed: none close.
    return f"{near} {100 * (1 - d / span):.1f}%", f"{nxt[0]} +{nxt[1]}" if nxt else "none close"


def call(tax, ties, d, span, rule, gap=None, cover=1.0):
    lab = merge_records([[tax[t][k] for k in RANKS] for t in ties])
    if d is None or span == 0:
        return "", "unassigned"
    ident = 1 - d / span
    allowed = None
    for rank, lim in RULES[rule]:
        if lim == "gap":
            ok = ident >= GAP_IDENTITY and gap is not None and gap >= GAP_EDITS
        else:
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
    tax = load_taxonomy(a.build)

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
        for m in cnt:
            cnt[m], f, d = denoise(cnt[m])
            stats[f"{m} folded"] += f
            stats[f"{m} noise"] += d
        per_run[run] = cnt
        for m in uniq:
            uniq[m].update(cnt[m])
    print(f"{stats['pairs']:,} pairs; {stats['primed']:,} with both primers; "
          f"{stats['merged']:,} merged; unique: "
          + ", ".join(f"{m} {len(c):,}" for m, c in uniq.items()), file=sys.stderr)
    print("denoising, reads folded into a commoner sequence / dropped as noise: "
          + ", ".join(f"{m} {stats[m + ' folded']:,} / {stats[m + ' noise']:,}" for m in uniq),
          file=sys.stderr)

    items = {m: list(c) for m, c in uniq.items()}
    jobs = []
    for m, its in items.items():
        # whole amplicon, or the two windows joined
        cand = minimap_candidates(os.path.join(a.build, "refdb.fasta"),
                                  its if m == "whole" else [x + "N" * 20 + y for x, y in its],
                                  a.workers)
        for i, s in enumerate(its):
            x, y = (s, "") if m == "whole" else s
            jobs.append(((m, i), "whole" if m == "whole" else "window", x, y, cand.get(i, [])))
    with Pool(a.workers, initializer=_init, initargs=(refs, tax)) as pool:
        res = {k: (d, sp, t, g, cv, nx) for k, d, sp, t, g, cv, nx in pool.imap_unordered(_score, jobs, chunksize=200)}

    idx = {m: {s: i for i, s in enumerate(its)} for m, its in items.items()}
    with open(a.prefix + ".reads.tsv", "w") as fh:
        fh.write("mode\tkey\treads\tedits\tspan\tidentity\tties\tnearest\tnext\n")
        for mode, its in items.items():
            for i, s in enumerate(its):
                d, sp, t, _, _, nx = res[(mode, i)]
                near, nxt = describe(tax, t, d, sp, nx)
                ident = "" if d is None or not sp else f"{1 - d / sp:.4f}"
                fh.write(f"{mode}\t{i}\t{uniq[mode][s]}\t{'' if d is None else d}\t{sp}\t{ident}\t{len(t)}\t"
                         f"{near}\t{nxt}\n")
    rows = collections.Counter()
    for run, cnt in per_run.items():
        for mode, c in cnt.items():
            for s, n in c.items():
                d, sp, t, g, cv, _ = res[(mode, idx[mode][s])]
                for rule in RULES:
                    rank, name = call(tax, t, d, sp, rule, g, cv)
                    rows[(run, mode, rule, rank or "unassigned", name)] += n
    with open(a.prefix + ".calls.tsv", "w") as fh:
        fh.write("run\tmode\trule\trank\tname\treads\n")
        for k, n in sorted(rows.items()):
            fh.write("\t".join(k) + f"\t{n}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
