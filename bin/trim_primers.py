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

  swapped  a pair with NC1 on read 2 and NC2 (or the motif) on read 1, as
          when both orientations were ligated (PRJNA1137257), is turned
          round first.
  read-through  a read longer than the amplicon runs on through the other
          primer: read 1 is cut where read 2 starts (the motif, read
          reverse-complemented), read 2 where read 1 starts (NC1, the same).

  clipped  a primer missing up to 3 bases at its 5' end (CGTCTGG... or
          GTCTGG... for NC1), as in a fifth of PRJNA1137257's pairs, is
          accepted.

Each of these acts only on a pair that needs it. --strict turns all
three off, which reproduces the trims of the real-data results before
2026-11-13 (they kept 0.4-0.6% fewer pairs).

--length R first truncates both reads to R bases, simulating a 2xR run
from longer reads (2x150 from 2x251).

--pretrimmed: the primers are already off, as submitted for PRJEB64582
(NC1 less its last two bases forward, ACGTCTGGTTAAGGGTTG, and a reverse
primer on the motif itself, ATGCTTAAGTTCAGCGGGTA): read 1 starts two
bases before the first ITS2 base, and those two are cut; read 2 starts
two bases past the motif and is kept whole. With no primer to go by, a
reversed pair (a third of PRJEB64582's, mostly whole runs deposited with
the two files swapped) is told by read 2's start, conserved 28S
(ATCACGACTGAGCTCA, <= 3 mismatches). --length R
counts the primers the reads would have carried (NC1, and NC2, its 7
bases and the 20-base primer), so 2x150 is as for the other sets.
"""
import argparse, collections, glob, gzip, multiprocessing as mp, os

NC1 = "ACGTCTGGTTCAGGGTTGTT"
NC2 = "TTAGTTTCTTTTCCTCCGCT"
MOTIF = "ATGCTTAAGTTCAGCGGG"
NC2_TAIL = 7            # 28S bases between NC2 and the motif
PRE2 = len(NC2) + NC2_TAIL + len(MOTIF) + 2   # read 2's start, pre-trimmed
R2_START = "ATCACGACTGAGCTCA"                     # ... and its first bases
RC = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s[::-1].translate(RC)


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


def primer_end(s, p, clip):
    """End of primer p near the read's start (<= 3 mismatches, 0-8 bases
    in), or, with clip, of p less up to 3 of its first bases; else None."""
    m = find(s, p, 3, 0, 8)
    if m:
        return m[0] + len(p)
    for k in range(1, 4 if clip else 1):
        m = find(s, p[k:], 3 - (k > 1), 0, 8)
        if m:
            return m[0] + len(p) - k
    return None


def read_through(s, p, mm, lo):
    """Where p first fits in s at or after lo, else len(s)."""
    m = find(s, p, mm, lo, len(s) - len(p))
    return m[0] if m else len(s)


def trim_run(job):
    raw1, out, length, pre, swap, through, clip = job
    run = os.path.basename(raw1).rsplit("_1.fastq.gz", 1)[0]
    raw2 = raw1.replace("_1.fastq.gz", "_2.fastq.gz")
    c = collections.Counter()
    with gzip.open(raw1, "rt") as a, gzip.open(raw2, "rt") as b, \
            gzip.open(os.path.join(out, f"{run}_1.fastq.gz"), "wt", compresslevel=4) as oa, \
            gzip.open(os.path.join(out, f"{run}_2.fastq.gz"), "wt", compresslevel=4) as ob:
        for (h1, s1, q1), (h2, s2, q2) in zip(records(a), records(b)):
            c["pairs"] += 1
            if pre:
                if swap and find(s1, R2_START, 3, 0, 3) and not find(s2, R2_START, 3, 0, 3):
                    h1, s1, q1, h2, s2, q2 = h2, s2, q2, h1, s1, q1
                    c["swapped"] += 1
                if length:
                    s1, q1 = s1[:length - len(NC1) + 2], q1[:length - len(NC1) + 2]
                    s2, q2 = s2[:length - PRE2], q2[:length - PRE2]
                cut1, cut2 = 2, 0
                end1 = read_through(s1, rc(MOTIF + "TA"), 2, 100) if through else len(s1)
                c["pre-trimmed"] += 1
            else:
                if length:
                    s1, q1, s2, q2 = s1[:length], q1[:length], s2[:length], q2[:length]
                if swap and not primer_end(s1, NC1, clip) and primer_end(s2, NC1, clip):
                    h1, s1, q1, h2, s2, q2 = h2, s2, q2, h1, s1, q1
                    c["swapped"] += 1
                cut1 = primer_end(s1, NC1, clip)
                if not cut1:
                    c["no NC1"] += 1
                    continue
                e = primer_end(s2, NC2, clip)
                if e:
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
                end1 = read_through(s1, rc(MOTIF), 2, cut1 + 100) if through else len(s1)
                end1 += len(MOTIF) if end1 < len(s1) else 0
            end2 = read_through(s2, rc(NC1), 3, cut2 + 100) if through else len(s2)
            c["read-through cut"] += end1 < len(s1) or end2 < len(s2)
            oa.write(f"{h1}{s1[cut1:end1]}\n+\n{q1[cut1:end1]}\n")
            ob.write(f"{h2}{s2[cut2:end2]}\n+\n{q2[cut2:end2]}\n")
            c["kept"] += 1
    return run, c


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("raw")
    ap.add_argument("out")
    ap.add_argument("--length", type=int, default=0)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--pretrimmed", action="store_true")
    ap.add_argument("--strict", action="store_true")
    # the defaults since 2026-11-13; accepted so older commands still run
    for flag in ("--swap", "--read-through", "--clipped"):
        ap.add_argument(flag, action="store_true", help=argparse.SUPPRESS)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    jobs = [(f, a.out, a.length, a.pretrimmed) + (not a.strict,) * 3 for f in sorted(glob.glob(os.path.join(a.raw, "*_1.fastq.gz")))]
    with mp.get_context("fork").Pool(a.workers) as pool:
        res = sorted(pool.imap_unordered(trim_run, jobs))
    cols = ["pairs", "kept", "read 2 at NC2", "read 2 in 28S", "no NC1", "no NC2 or 28S",
            "swapped", "pre-trimmed", "read-through cut"]
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
