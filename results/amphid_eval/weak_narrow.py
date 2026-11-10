"""Species share if 'other species' counted only full-length or multi-accession records.
weak.py <acc dir> <n top sequences>"""
import csv, glob, sys, collections
from multiprocessing import Pool
sys.path.insert(0, "/Users/sd21/Documents/WORK/informatics_workbook/nf-amphid-pipeline/bin")
import assign_db as adb, ampseqlib as al
B = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/work/build_v0.1.0"
refs = dict(al.read_fasta(B + "/refdb.fasta")); tax = al.load_taxonomy(B + "/taxonomy.tsv")
nmem = collections.Counter(r["id"] for r in csv.DictReader(open(B + "/members.tsv"), delimiter="\t") if r["role"] == "kept")
FULL = ("complete", "near_complete")
import json, os
weak = {i for i, t in tax.items() if t.get("grade") not in FULL and nmem[i] == 1}
SETS = {k: set(v) for k, v in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "narrow_sets.json"))).items()}
SETS["all weak"] = weak
def calls(s):
    x, _, y = s.partition("N" * 10)
    parts = [(x, "start"), (y, "end")] if y else [(x, "whole")]
    sc = adb.span_scores(parts, list(refs), refs, 47)
    a = adb.span_call(sc, tax, "gap", adb.GAP_FLOORS)
    # weak records may still be the best hit; they just don't count as the nearest other species
    best_species = {tax[c]["species"] for c in a["ties"]}
    out = {k: adb.span_call([t for t in sc if t[2] not in v or tax[t[2]]["species"] in best_species],
                            tax, "gap", adb.GAP_FLOORS)["rank"] for k, v in SETS.items()}
    return s, a["rank"], out
if __name__ == "__main__":
    c = collections.Counter()
    for p in glob.glob(sys.argv[1] + "/*.clusters.tsv"):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["origin"] != "chimera": c[r["consensus"]] += int(r["reads"])
    top = [s for s, _ in c.most_common(int(sys.argv[2]))]
    with Pool(8) as p: res = p.map(calls, top)
    tot = sum(c.values()); cov = sum(c[s] for s in top)
    a = sum(c[s] for s, ra, _o in res if ra == "species")
    print(f"{sys.argv[1].split('/')[-1]}: top {len(top)} = {100*cov/tot:.1f}% of reads; species as released {100*a/tot:.1f}%; "
          + "; ".join(f"{k} {100*sum(c[s] for s, _r, o in res if o[k] == 'species')/tot:.1f}%" for k in SETS))
