#!/usr/bin/env python3
"""Measure, trim and grade the ITS2 records bin/its2.py could not anchor.

    bin/its2_partials.py <its2 cut prefix> <out prefix> [threads]

Two thirds of the backbone's source records lack NC1, the 28S stretch or
both: submitters trimmed at a primer, at the ITS2 boundary, or somewhere
inside it. Each such record is mapped against the records anchored at
both ends ("complete"), and its best relative stands in for the missing
anchor:

  trim    bases running past the relative's ITS2 are cut off - 5.8S
          before NC1 (not in the reads at all) or a stub of the 28S
          stretch too short for the anchor to find. Only where the
          relative's own end is reached (within END_SLOP bases) and the
          overhang is at least MIN_OVERHANG, so a divergent tail that
          merely fails to align is left alone.
  grade   coverage of the relative's ITS2 after trimming:
            complete        both ends anchored by bin/its2.py
            near_complete   >= 98%
            partial_90      90-98%
            partial_75      75-90%
            partial         < 75%
            unmeasured      no complete relative at >= 80% identity

80% is low on purpose: sister species within a genus are often 85-90%
identical over ITS2 and still align end to end, which is all a boundary
estimate needs.

minimap2 runs with -f 1000000: its default drops minimizers shared by
too many targets, and thousands of near-identical strongylid records
made it miss whole species (Strongylus vulgaris).

Writes <prefix>.fasta (every record, trimmed) and <prefix>.tsv
(accession, grade, relative, identity, coverage, missing/trimmed bases
at each end, length).
"""

from __future__ import annotations

import collections
import csv
import os
import subprocess
import sys
import tempfile

MIN_IDENTITY = 0.80
MIN_OVERHANG = 8
END_SLOP = 3
ANCHORED = ("primer", "loose", "part")


def read_fasta(path):
    seqs, cur = {}, None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def grade(cov):
    return ("near_complete" if cov >= 0.98 else "partial_90" if cov >= 0.90
            else "partial_75" if cov >= 0.75 else "partial")


def main() -> int:
    cut, prefix = sys.argv[1:3]
    threads = sys.argv[3] if len(sys.argv) > 3 else "8"
    meta = {r["accession"]: r for r in csv.DictReader(open(cut + ".tsv"), delimiter="\t")}
    seqs = read_fasta(cut + ".fasta")
    complete = {a for a, r in meta.items() if r["start"] in ANCHORED and r["end"] in ANCHORED}

    with tempfile.TemporaryDirectory() as d:
        cf, pf = os.path.join(d, "c.fa"), os.path.join(d, "p.fa")
        with open(cf, "w") as c, open(pf, "w") as p:
            for a, s in seqs.items():
                (c if a in complete else p).write(f">{a}\n{s}\n")
        paf = subprocess.run(["minimap2", "-c", "-k", "11", "-w", "5", "-f", "1000000",
                              "-N", "50", "-p", "0.5", "-t", threads, "--secondary=yes", cf, pf],
                             capture_output=True, text=True, check=True).stdout
    best = {}
    for line in paf.splitlines():
        f = line.split("\t")
        if f[4] != "+":
            continue
        m = int(f[9])
        if f[0] not in best or m > best[f[0]][0]:
            best[f[0]] = (m, int(f[10]), int(f[1]), int(f[2]), int(f[3]),
                          f[5], int(f[6]), int(f[7]), int(f[8]))

    rows, tally = [], collections.Counter()
    with open(prefix + ".fasta", "w") as fa:
        for a, s in seqs.items():
            row = {"accession": a, "grade": "complete", "relative": "", "identity": "",
                   "coverage": "", "missing_start": "", "missing_end": "",
                   "trimmed_start": 0, "trimmed_end": 0}
            if a not in complete:
                b = best.get(a)
                if not b or b[0] / b[1] < MIN_IDENTITY:
                    row["grade"] = "unmeasured"
                else:
                    m, al, ql, qs, qe, rel, tl, ts, te = b
                    if ts <= END_SLOP and qs >= MIN_OVERHANG:
                        row["trimmed_start"] = qs
                    if tl - te <= END_SLOP and ql - qe >= MIN_OVERHANG:
                        row["trimmed_end"] = ql - qe
                    s = s[row["trimmed_start"]:len(s) - row["trimmed_end"]]
                    cov = (te - ts) / tl
                    row.update(grade=grade(cov), relative=rel, identity=f"{m / al:.3f}",
                               coverage=f"{cov:.3f}", missing_start=ts, missing_end=tl - te)
            row["length"] = len(s)
            fa.write(f">{a}\n{s}\n")
            rows.append(row)
            tally[row["grade"]] += 1
            tally["trimmed"] += bool(row["trimmed_start"] or row["trimmed_end"])
    with open(prefix + ".tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for k, n in tally.most_common():
        print(f"{n:6d}  {k}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
