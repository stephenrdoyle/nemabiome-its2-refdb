#!/usr/bin/env python3
"""Primer-trim NC1/NC2 read pairs for DADA2, both read 2 layouts.

    bin/trim_primers.py <raw dir> <out dir> [--length R] [--workers N]

<raw dir> holds <run>_1.fastq.gz / <run>_2.fastq.gz. Writes the trimmed
pairs to <out dir> under the same names, and <out dir>/trim_stats.tsv.

  read 1  NC1 (ACGTCTGGTTCAGGGTTGTT), <= 3 mismatches, starting within the
          first 8 bases (0-3 N spacers and a little slack); cut after it.
          No NC1: the pair is dropped.
  read 2  starts either with NC2 (TTAGTTTCTTTTCCTCCGCT, same tolerance),
          then 7 bases of 28S, or directly inside the 28S stretch, as in
          half the runs of PRJNA1524548, which a trim requiring NC2 throws
          away. Either way it is cut at the conserved 28S motif
          ATGCTTAAGTTCAGCGGG (<= 2 mismatches, within 12 bases of NC2's end
          or of the read's start), so every read 2 starts on the same base;
          the 7 bases dropped are 28S, which no record holds. NC2 without
          the motif: cut 7 bases past NC2. Neither: the pair is dropped.

--length R first truncates both reads to R bases, simulating a 2xR run
from longer reads (2x150 from 2x251).
"""
import argparse, collections, glob, gzip, multiprocessing as mp, os

NC1 = "ACGTCTGGTTCAGGGTTGTT"
NC2 = "TTAGTTTCTTTTCCTCCGCT"
MOTIF = "ATGCTTAAGTTCAGCGGG"
NC2_TAIL = 7            # 28S bases between NC2 and the motif


def find(s, p, mm, lo, hi):
    """(offset, mismatches) of the best Hamming fit of p starting in [lo, hi]."""
    best = None
    for off in range(lo, hi + 1):
        w = s[off:off + len(p)]
        if len(w) < len(p):
            break
        d = sum(a != b for a, b in zip(w, p))
        if d <= mm and (best is None or d < best[1]):
            best = (off, d)
    return best


def records(fh):
    while True:
        h = fh.readline()
        if not h:
            return
        s, _, q = fh.readline().rstrip("\n"), fh.readline(), fh.readline().rstrip("\n")
        yield h, s, q


def trim_run(job):
    raw1, out, length = job
    run = os.path.basename(raw1).rsplit("_1.fastq.gz", 1)[0]
    raw2 = raw1.replace("_1.fastq.gz", "_2.fastq.gz")
    c = collections.Counter()
    with gzip.open(raw1, "rt") as a, gzip.open(raw2, "rt") as b, \
            gzip.open(os.path.join(out, f"{run}_1.fastq.gz"), "wt", compresslevel=4) as oa, \
            gzip.open(os.path.join(out, f"{run}_2.fastq.gz"), "wt", compresslevel=4) as ob:
        for (h1, s1, q1), (h2, s2, q2) in zip(records(a), records(b)):
            c["pairs"] += 1
            if length:
                s1, q1, s2, q2 = s1[:length], q1[:length], s2[:length], q2[:length]
            p1 = find(s1, NC1, 3, 0, 8)
            if not p1:
                c["no NC1"] += 1
                continue
            p2 = find(s2, NC2, 3, 0, 8)
            if p2:
                e = p2[0] + len(NC2)
                m = find(s2, MOTIF, 2, e, e + 12)
                cut2 = m[0] if m else e + NC2_TAIL
                c["read 2 at NC2"] += 1
            else:
                m = find(s2, MOTIF, 2, 0, 12)
                if not m:
                    c["no NC2 or 28S"] += 1
                    continue
                cut2 = m[0]
                c["read 2 in 28S"] += 1
            cut1 = p1[0] + len(NC1)
            oa.write(f"{h1}{s1[cut1:]}\n+\n{q1[cut1:]}\n")
            ob.write(f"{h2}{s2[cut2:]}\n+\n{q2[cut2:]}\n")
            c["kept"] += 1
    return run, c


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("raw")
    ap.add_argument("out")
    ap.add_argument("--length", type=int, default=0)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = [(f, a.out, a.length) for f in sorted(glob.glob(os.path.join(a.raw, "*_1.fastq.gz")))]
    with mp.get_context("fork").Pool(a.workers) as pool:
        res = sorted(pool.imap_unordered(trim_run, jobs))
    cols = ["pairs", "kept", "read 2 at NC2", "read 2 in 28S", "no NC1", "no NC2 or 28S"]
    tot = collections.Counter()
    with open(os.path.join(a.out, "trim_stats.tsv"), "w") as fh:
        fh.write("run\t" + "\t".join(c.replace(" ", "_") for c in cols) + "\n")
        for run, c in res:
            tot.update(c)
            fh.write(run + "\t" + "\t".join(str(c[x]) for x in cols) + "\n")
    print(f"{a.raw} -> {a.out}: {len(res)} runs; " + ", ".join(
        f"{x} {tot[x]:,}" + (f" ({100 * tot[x] / tot['pairs']:.1f}%)" if x != "pairs" else "") for x in cols))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
