"""score.py <acc dir>: read-weighted species share, pooled; plus fates."""
import csv, glob, json, sys, collections
d = sys.argv[1]
tot = sp = ge = 0
fate = collections.Counter(); comp = collections.Counter()
for p in glob.glob(f"{d}/*.clusters.tsv"):
    for r in csv.DictReader(open(p), delimiter="\t"):
        n = int(r["reads"]); tot += n
        if r["origin"] == "chimera": comp["(chimera)"] += n; continue
        sp += n * (r["rank"] == "species"); ge += n * (r["rank"] == "genus")
        comp[r["taxon"] if r["rank"] == "species" else f"({r['rank']}) {r['taxon']}"] += n
for p in glob.glob(f"{d}/*.assign_db.json"):
    fate.update(json.load(open(p))["reads"])
print(f"{d.split('/')[-1]}: {len(glob.glob(d+'/*.clusters.tsv'))} samples, {tot:,} clustered pairs; species {100*sp/tot:.1f}%  genus {100*ge/tot:.1f}%  other {100*(tot-sp-ge)/tot:.1f}%")
print("  fates:", dict(fate))
print("  top:", ", ".join(f"{k} {100*v/tot:.1f}%" for k, v in comp.most_common(9)))
