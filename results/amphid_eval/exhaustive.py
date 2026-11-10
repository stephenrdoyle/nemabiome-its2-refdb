import csv, glob, sys, collections, os
from multiprocessing import Pool
sys.path.insert(0, "/Users/sd21/Documents/WORK/informatics_workbook/nf-amphid-pipeline/bin")
import assign_db as adb, ampseqlib as al
R = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/work/build_v0.1.0"
refs = dict(al.read_fasta(R + "/refdb.fasta"))
tax = al.load_taxonomy(R + "/taxonomy.tsv")
def one(s):
    x, _, y = s.partition("N" * 10)
    parts = [(x, "start"), (y, "end")] if y else [(x, "whole")]
    c = adb.span_call(adb.span_scores(parts, list(refs), refs, 47), tax, "gap", adb.GAP_FLOORS)
    return s, (c["rank"], c["gap"])
if __name__ == "__main__":
    c = collections.Counter(); info = {}
    for p in glob.glob(sys.argv[1] + "/*.clusters.tsv"):
        for r in csv.DictReader(open(p), delimiter="\t"):
            if r["origin"] != "chimera":
                c[r["consensus"]] += int(r["reads"]); info[r["consensus"]] = (r["rank"], r["gap"])
    top = [s for s, _ in c.most_common(int(sys.argv[2]))]
    with Pool(8) as pool:
        res = dict(pool.map(one, top))
    tot = sum(c[s] for s in top)
    same = sum(c[s] for s in top if res[s][0] == info[s][0])
    print(f"top {len(top)} sequences, {tot:,} reads ({100*tot/sum(c.values()):.1f}%): same rank exhaustive vs amphid on {100*same/tot:.2f}% of reads")
    for s in top:
        if res[s][0] != info[s][0]: print(c[s], "amphid", info[s], "exhaustive", res[s])
