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
weak = {i for i, t in tax.items() if t.get("grade") not in FULL and nmem[i] == 1}
def calls(s):
    x, _, y = s.partition("N" * 10)
    parts = [(x, "start"), (y, "end")] if y else [(x, "whole")]
    sc = adb.span_scores(parts, list(refs), refs, 47)
    a = adb.span_call(sc, tax, "gap", adb.GAP_FLOORS)
    # weak records may still be the best hit; they just don't count as the nearest other species
    best_species = {tax[c]["species"] for c in a["ties"]}
    sc2 = [t for t in sc if t[2] not in weak or tax[t[2]]["species"] in best_species]
    b = adb.span_call(sc2, tax, "gap", adb.GAP_FLOORS)
    return s, a["rank"], b["rank"], a["next"]
if __name__ == "__main__":
    c = collections.Counter()
    for p in glob.glob(sys.argv[1] + "/*.clusters.tsv"):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["origin"] != "chimera": c[r["consensus"]] += int(r["reads"])
    top = [s for s, _ in c.most_common(int(sys.argv[2]))]
    with Pool(8) as p: res = p.map(calls, top)
    tot = sum(c.values()); cov = sum(c[s] for s in top)
    a = sum(c[s] for s, ra, rb, _ in res if ra == "species"); b = sum(c[s] for s, ra, rb, _ in res if rb == "species")
    print(f"{sys.argv[1].split('/')[-1]}: top {len(top)} sequences = {100*cov/tot:.1f}% of reads; species {100*a/tot:.1f}% as released -> {100*b/tot:.1f}% if single-accession partials are not 'another species' ({len(weak):,} such records in the release)")
