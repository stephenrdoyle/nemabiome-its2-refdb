#!/usr/bin/env python3
"""Would a species rule set per genus or family beat one rule for all?

    bin/loo_group_rules.py <build dir> <leave-one-out .tsv> > <out.txt>

From bin/its2_leave_one_out.py's raw scores; nothing is rescored.

A species call needs one best species, identity >= floor and gap >= g.
Global: one (floor, g) for all. Per group: each genus (or family) of the
best species gets its own (floor, g), chosen on the other folds' queries
to maximise right - lam * wrong (present wrong + hidden called); groups
with < MIN_N training queries fall back to family, then to global.
Species split into 5 folds, so a species' queries are never in the
training set used for them. lam is the cost of one wrong call in right
calls; the current rule (99%, gap 4) is printed first for comparison.
"""
import collections, csv, random, sys

FLOORS = [0.99, 0.985, 0.98, 0.97]
GAPS = [1, 2, 3, 4, 6, 10]
MIN_N = 8

build, loo = sys.argv[1:3]
fam = {}
for r in csv.DictReader(open(f"{build}/taxonomy.tsv"), delimiter="\t"):
    if r["species"]:
        fam[r["species"]] = (r["genus"], r["family"])
rows = []
for r in csv.DictReader(open(loo), delimiter="\t"):
    ties = r["tied_species"].split(";") if r["tied_species"] else []
    ok = r["identity"] != "" and r["gap"] != "" and len(ties) == 1
    rows.append(dict(sp=r["species"], mode=r["mode"], view=r["view"], ok=ok,
                     ident=float(r["identity"]) if ok else 0, gap=int(r["gap"]) if ok else 0,
                     best=ties[0] if ok else "", base=r["gap_rank"] == "species",
                     grp=fam.get(ties[0], ("", "")) if ok else ("", "")))


def outcome(q, rule):
    f, g = rule
    if not (q["ok"] and q["ident"] >= f and q["gap"] >= g):
        return 0, 0
    right = q["mode"] == "present" and q["best"] == q["sp"]
    return (1, 0) if right else (0, 1)


def best_rule(qs, lam):
    score = {}
    for f in FLOORS:
        for g in GAPS:
            r = w = 0
            for q in qs:
                a, b = outcome(q, (f, g))
                r += a; w += b
            score[f, g] = r - lam * w
    return max(score, key=lambda k: (score[k], k[0], k[1]))


def evaluate(rule_of, view):
    c = collections.Counter()
    for q in rows:
        if q["view"] != view:
            continue
        a, b = outcome(q, rule_of(q))
        c[q["mode"], "n"] += 1
        c[q["mode"], "right"] += a
        c[q["mode"], "wrong"] += b
    return c


species = sorted({q["sp"] for q in rows})
random.Random(3).shuffle(species)
fold = {s: i % 5 for i, s in enumerate(species)}


def line(lab, c):
    pn, hn = c["present", "n"], c["hidden", "n"]
    return (f"  {lab:34s} present right {100 * c['present', 'right'] / pn:5.1f}%  "
            f"wrong {100 * c['present', 'wrong'] / pn:4.1f}%   hidden wrong {100 * c['hidden', 'wrong'] / hn:4.1f}%")


for view in ("whole", "2x150"):
    print(f"== {view}")
    print(line("current (99%, gap 4)", evaluate(lambda q: (0.99, 4), view)))
    for lam in (2, 4, 8):
        res = {}
        for level in ("global", "genus", "family", "genus>family"):
            chosen = {}
            for k in range(5):
                train = [q for q in rows if q["view"] == view and fold[q["sp"]] != k]
                glob = best_rule(train, lam)
                by = {0: collections.defaultdict(list), 1: collections.defaultdict(list)}
                for q in train:
                    if q["ok"]:
                        by[0][q["grp"][0]].append(q); by[1][q["grp"][1]].append(q)
                rules = {0: {}, 1: {}}
                for i in (0, 1):
                    for key, qs in by[i].items():
                        if len(qs) >= MIN_N:
                            rules[i][key] = best_rule(qs, lam)
                chosen[k] = (glob, rules)

            def rule_of(q, level=level):
                glob, rules = chosen[fold[q["sp"]]]
                if level == "global":
                    return glob
                if level in ("genus", "genus>family") and q["grp"][0] in rules[0]:
                    return rules[0][q["grp"][0]]
                if level in ("family", "genus>family") and q["grp"][1] in rules[1]:
                    return rules[1][q["grp"][1]]
                return glob
            print(line(f"lam {lam}: {level}", evaluate(rule_of, view)))
    print()
