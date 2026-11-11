#!/usr/bin/env python3
"""Is an exact match with a sister species close by safe to call?

    bin/loo_gap_split.py <leave-one-out .tsv> > <out.txt>

From bin/its2_leave_one_out.py's raw scores; nothing is rescored. Queries
whose best hit is one species (no tie), by identity band and gap (edits
from the best species to the next):

  present  the query's own record is removed, its species stays: a
           species call is right if the best species is its own
  hidden   its whole species is removed: any species call is wrong, the
           over-classification the gap rule guards against

Then the rule with a smaller gap for exact matches only (identity 100%,
the gap still 4 below it), as recall, precision and over-classification
over all queries, as in bin/loo_metrics.py.
"""
import collections, csv, sys

BANDS = [("100%", 1.0, 1.01), ("99-100%", 0.99, 1.0)]
GAPS = [(1, 1, "1"), (2, 2, "2"), (3, 3, "3"), (4, 4, "4"), (5, 6, "5-6"), (7, 10 ** 6, "7+")]

rows = []
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    ties = r["tied_species"].split(";") if r["tied_species"] else []
    ok = r["identity"] != "" and r["gap"] != "" and len(ties) == 1
    rows.append(dict(sp=r["species"], mode=r["mode"], view=r["view"], ok=ok,
                     ident=float(r["identity"]) if ok else 0, gap=int(r["gap"]) if ok else 0,
                     best=ties[0] if ok else ""))


def named(q, exact_gap):
    if not q["ok"] or q["ident"] < 0.99:
        return False
    return q["gap"] >= (exact_gap if q["ident"] >= 1.0 else 4)


for view in ("whole", "2x150"):
    vq = [q for q in rows if q["view"] == view]
    npres = sum(q["mode"] == "present" for q in vq)
    nhid = sum(q["mode"] == "hidden" for q in vq)
    print(f"== {view}: {npres:,} present, {nhid:,} hidden queries")
    print(f"   {'identity':9s} {'gap':>4s} {'present':>8s} {'right':>7s} {'wrong':>6s} {'hidden':>7s}"
          f"  {'right if called':>15s}")
    for blab, lo, hi in BANDS:
        for g0, g1, glab in GAPS:
            cell = [q for q in vq if q["ok"] and lo <= q["ident"] < hi and g0 <= q["gap"] <= g1]
            p = [q for q in cell if q["mode"] == "present"]
            h = [q for q in cell if q["mode"] == "hidden"]
            right = sum(q["best"] == q["sp"] for q in p)
            share = f"{100 * right / len(p):5.1f}% of present" if p else ""
            print(f"   {blab:9s} {glab:>4s} {len(p):8,d} {right:7,d} {len(p) - right:6,d} {len(h):7,d}  {share}")
    print()
    print("   Exact matches named at a smaller gap (99-100%: still gap 4):")
    for eg in (4, 3, 2, 1):
        c = collections.Counter()
        for q in vq:
            if named(q, eg):
                c[q["mode"]] += 1
                c["right"] += q["mode"] == "present" and q["best"] == q["sp"]
        lab = "current" if eg == 4 else f"exact at gap >= {eg}"
        print(f"   {lab:22s} recall {100 * c['right'] / npres:5.1f}%  precision "
              f"{100 * c['right'] / c['present']:5.1f}%  over-class. {100 * c['hidden'] / nhid:5.1f}%")
    print()
