#!/usr/bin/env python3
"""Does a distance-scaled species rule below 99% hold up on held-out species?

    bin/loo_scaled_floor.py <leave-one-out .tsv> > <out.txt>

From bin/its2_leave_one_out.py's raw scores; nothing is rescored.

The current rule names a species at >= 99% identity and gap >= 4. The
scaled rule keeps that, and below 99% also names it at identity >= F
when the gap is >= max(b, m x edits to the best): the further the read is
from its best species, the further the next species must be. "None"
(no branch below 99%, the current rule) is in the grid. Compared
with one plain (floor, gap) rule. Each rule's parameters are chosen on
four fifths of the species to maximise right - lam * wrong (present
wrong + hidden called) and scored on the fifth held out, as in
bin/loo_group_rules.py (same folds).

Recall = present right / present queries; precision = present right /
present called; over-classification = hidden called / hidden queries.
"""
import collections, csv, itertools, random, sys

FLOORS = [0.99, 0.985, 0.98, 0.97]
GAPS = [1, 2, 3, 4, 6, 10]
S_FLOORS = [0.985, 0.98, 0.975, 0.97, 0.96, 0.95]
S_BASE = [4, 6, 8, 10]
S_MULT = [0, 1, 1.5, 2, 3]

rows = []
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    ties = r["tied_species"].split(";") if r["tied_species"] else []
    ok = r["identity"] != "" and r["gap"] != "" and len(ties) == 1
    rows.append(dict(sp=r["species"], mode=r["mode"], view=r["view"], ok=ok,
                     ident=float(r["identity"]) if ok else 0, gap=int(r["gap"]) if ok else 0,
                     edits=int(r["edits"]) if ok else 0, best=ties[0] if ok else ""))


def plain(q, rule):
    f, g = rule
    return q["ok"] and q["ident"] >= f and q["gap"] >= g


def scaled(q, rule):
    f, b, m = rule
    if not q["ok"]:
        return False
    if q["ident"] >= 0.99 and q["gap"] >= 4:
        return True
    return f is not None and q["ident"] >= f and q["gap"] >= max(b, m * q["edits"])


def tally(qs, named, rule):
    c = collections.Counter()
    for q in qs:
        c[q["mode"], "n"] += 1
        if named(q, rule):
            c[q["mode"], "called"] += 1
            c[q["mode"], "right"] += q["mode"] == "present" and q["best"] == q["sp"]
    return c


def best(qs, named, grid, lam):
    def score(rule):
        c = tally(qs, named, rule)
        wrong = c["present", "called"] - c["present", "right"] + c["hidden", "called"]
        return c["present", "right"] - lam * wrong
    return max(grid, key=lambda r: (score(r), r))


species = sorted({q["sp"] for q in rows})
random.Random(3).shuffle(species)
fold = {s: i % 5 for i, s in enumerate(species)}


def line(lab, c):
    pn, hn, pc = c["present", "n"], c["hidden", "n"], c["present", "called"]
    return (f"  {lab:30s} recall {100 * c['present', 'right'] / pn:5.1f}%  "
            f"precision {100 * c['present', 'right'] / pc:5.1f}%  "
            f"over-class. {100 * c['hidden', 'called'] / hn:5.1f}%")


plain_grid = list(itertools.product(FLOORS, GAPS))
scaled_grid = [(None, 0, 0)] + list(itertools.product(S_FLOORS, S_BASE, S_MULT))
for view in ("whole", "2x150"):
    vq = [q for q in rows if q["view"] == view]
    print(f"\n== {view}")
    print(line("current (99%, gap 4)", tally(vq, plain, (0.99, 4))))
    for lam in (2, 4, 8):
        for lab, named, grid in (("plain", plain, plain_grid), ("scaled", scaled, scaled_grid)):
            total, picks = collections.Counter(), []
            for k in range(5):
                rule = best([q for q in vq if fold[q["sp"]] != k], named, grid, lam)
                picks.append(rule)
                total.update(tally([q for q in vq if fold[q["sp"]] == k], named, rule))
            print(line(f"lam {lam}: {lab}", total) + "   chosen " +
                  ", ".join(sorted({"/".join(str(x) for x in p) for p in picks})))
