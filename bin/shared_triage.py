#!/usr/bin/env python3
"""Which name does a shared ITS2 belong to?

    bin/shared_triage.py <build dir> > <out.tsv>

Where one ITS2 carries two species' names - a shared record (members
naming two or more species) or two named records 0 edits apart over
their span (species_margins.tsv) - ask, for each name, how the sequence
sits among that species' other records:

  n        other records of the species
  near     edits to the nearest of them
  within3  share of them within 3 edits
  spread   the species' own median nearest-neighbour distance, the yardstick
           for "typical" (a species with diverse ITS2 forgives more)

and how far up the hierarchy the two names part (species = same genus,
genus = same family, family, order). Then a call:

  fits A       typical of A (near <= max(3, 2 x spread)), not of B, and B
               has other records to judge by: the B label is the likely
               mislabel
  both         typical of both: a real shared ITS2 (sister species, or
               morphs under two names)
  neither      typical of neither: an odd sequence, check its source
  one-sided    a name with no other records: nothing to judge by; the
               higher the split, the more suspect

Distances are bin/build.py's span_distance (a partial compared over the
bases it has). Members are counted by accession per name: two
independent sources on one side outweigh one.
"""
import collections, csv, statistics, sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from build import span_distance  # noqa: E402

RANKS = ["order", "family", "genus", "species"]

b = sys.argv[1]
seq, k = {}, None
for line in open(f"{b}/refdb.fasta"):
    line = line.strip()
    if line.startswith(">"):
        k = line[1:].split()[0]
        seq[k] = []
    else:
        seq[k].append(line)
seq = {k: "".join(v) for k, v in seq.items()}
tax = {r["id"]: r for r in csv.DictReader(open(f"{b}/taxonomy.tsv"), delimiter="\t")}
mem = collections.defaultdict(list)
for r in csv.DictReader(open(f"{b}/members.tsv"), delimiter="\t"):
    if r["id"] and r["role"] != "dropped":
        mem[r["id"]].append(r)
margins = {r["id"]: r for r in csv.DictReader(open(f"{b}/species_margins.tsv"), delimiter="\t")}

lineage = {}                                  # species -> (order, family, genus)
by_sp = collections.defaultdict(list)         # species -> its named records
for i, t in tax.items():
    if t["species"] and t["species_status"] == "named":
        by_sp[t["species"]].append(i)
        lineage[t["species"]] = (t["order"], t["family"], t["genus"])
for rs in mem.values():
    for r in rs:
        p = r["label"].split(";")
        if len(p) >= 8 and p[7] and " x " not in p[7]:
            lineage.setdefault(p[7], (p[4], p[5], p[6]))

_spread = {}


def spread(sp):
    if sp not in _spread:
        ids = by_sp.get(sp, [])[:40]
        nn = [min(span_distance(seq[i], seq[j]) for j in ids if j != i) for i in ids] if len(ids) > 1 else []
        _spread[sp] = statistics.median(nn) if nn else None
    return _spread[sp]


def profile(s, sp, exclude):
    others = [j for j in by_sp.get(sp, []) if j not in exclude]
    if not others:
        return dict(n=0, near="", within3="", spread="", typical=None)
    d = [span_distance(s, seq[j]) for j in others]
    sp_ = spread(sp)
    near = min(d)
    return dict(n=len(others), near=near, within3=f"{sum(x <= 3 for x in d) / len(d):.2f}",
                spread="" if sp_ is None else f"{sp_:g}",
                typical=near <= max(3, 2 * (sp_ or 0)))


def split_rank(a, b):
    la, lb = lineage.get(a, ("", "", "")), lineage.get(b, ("", "", ""))
    for r, x, y in zip(RANKS[:3], la, lb):
        if x != y:
            return r
    return "species"


cases = []
for i, t in tax.items():                       # shared records
    if t["species_status"] != "shared":
        continue
    names = collections.Counter()
    for r in mem[i]:
        p = r["label"].split(";")
        if len(p) >= 8 and p[7] and " x " not in p[7] and " sp." not in p[7]:
            names[p[7]] += 1
    if len(names) >= 2:
        (a, na), (bb, nb) = names.most_common(2)
        cases.append((i, i, a, bb, na, nb))
seen = set()
for i, r in margins.items():                   # named records 0 edits apart
    n = r["nearest_whole"]
    if r["margin_whole"] == "0" and n and frozenset((i, n)) not in seen:
        seen.add(frozenset((i, n)))
        if tax[i]["species"] and tax[n]["species"] and tax[i]["species"] != tax[n]["species"]:
            cases.append((i, n, tax[i]["species"], tax[n]["species"], len(mem[i]), len(mem[n])))

cols = ["record", "other_record", "A", "B", "split", "acc_A", "acc_B",
        "n_A", "near_A", "within3_A", "spread_A", "n_B", "near_B", "within3_B", "spread_B",
        "call", "accessions_A", "accessions_B"]
print("\t".join(cols))
for i, j, a, bb, na, nb in cases:
    s = seq[i]
    pa, pb = profile(s, a, {i, j}), profile(s, bb, {i, j})
    if pa["typical"] is None or pb["typical"] is None:
        call = "one-sided"
    elif pa["typical"] and pb["typical"]:
        call = "both"
    elif pa["typical"]:
        call = "fits A"
    elif pb["typical"]:
        call = "fits B"
    else:
        call = "neither"
    acc = lambda rid, sp: ",".join(r["accession"] for r in mem[rid] if r["label"].endswith(";" + sp))[:80]
    print("\t".join(map(str, [i, j if j != i else "", a, bb, split_rank(a, bb), na, nb,
                              pa["n"], pa["near"], pa["within3"], pa["spread"],
                              pb["n"], pb["near"], pb["within3"], pb["spread"], call,
                              acc(i, a) or acc(j, a), acc(j, bb) or acc(i, bb)])))
