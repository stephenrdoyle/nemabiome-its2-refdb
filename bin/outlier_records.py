#!/usr/bin/env python3
"""Single records that crowd another species.

    bin/outlier_records.py <build dir> > <out.txt>

From the build's species_margins.tsv (nearest other-species record and
its distance, read whole). For each record R, the records of another
species X whose nearest other-species record is R, within 3 edits: R
holds them under the gap rule. Listed when R does this for >= 3 records
of X, with how many other records of R's own species are also within 3
edits of X ("sib"). sib = 0 means R is the odd one out of its species:
a candidate mislabel, or a partial that ends where X differs. Unnamed
records (hybrids, cf.) are listed too, though the classifier ignores
them.
"""
import collections, csv, sys

b = sys.argv[1]
tax = {r["id"]: r for r in csv.DictReader(open(f"{b}/taxonomy.tsv"), delimiter="\t")}
mem = collections.defaultdict(list)
for r in csv.DictReader(open(f"{b}/members.tsv"), delimiter="\t"):
    if r["id"] and r["role"] != "dropped":
        mem[r["id"]].append(r["accession"])
mg = {r["id"]: r for r in csv.DictReader(open(f"{b}/species_margins.tsv"), delimiter="\t")
      if r["margin_whole"] != ""}
sp = lambda i: tax[i]["species"] or f"({tax[i]['genus']} {tax[i]['species_status']})"
held = collections.defaultdict(list)
for i, r in mg.items():
    if int(r["margin_whole"]) <= 3 and sp(r["nearest_whole"]) != sp(i):
        held[r["nearest_whole"]].append(i)
nsp = collections.Counter(sp(i) for i in mg)
rows = []
for R, lst in held.items():
    X = collections.Counter(sp(i) for i in lst).most_common(1)[0][0]
    n = sum(sp(i) == X for i in lst)
    if n < 3:
        continue
    sib = sum(1 for j in mg if sp(j) == sp(R) and j != R and sp(mg[j]["nearest_whole"]) == X
              and int(mg[j]["margin_whole"]) <= 3)
    rows.append((n, nsp[X], R, sp(R), nsp[sp(R)], sib, X, tax[R]["grade"], ",".join(mem[R][:3])))
print("held X / X records | record R | R's species (records) | sib | X | grade | accessions")
for r in sorted(rows, reverse=True):
    print(f"{r[0]:3d}/{r[1]:<4d} {r[2]} {r[3]} ({r[4]}) sib={r[5]} -> {r[6]}  {r[7]} [{r[8]}]")
