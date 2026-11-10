import csv, collections, sys
from multiprocessing import Pool
sys.path.insert(0, "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/bin")
from build import span_distance, read_fasta
B = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/work/build_v0.1.0"
seqs = read_fasta(B + "/refdb.fasta")
tax = {r["id"]: r for r in csv.DictReader(open(B + "/taxonomy.tsv"), delimiter="\t")}
by = collections.defaultdict(list)
for i, t in tax.items():
    if t["species"]: by[t["species"]].append(i)
def nearest_same(args):
    i, ids = args
    best = None
    for j in ids:
        if j != i:
            d = span_distance(seqs[i], seqs[j])
            if best is None or d < best: best = d
            if best == 0: break
    return i, best
if __name__ == "__main__":
    n = sorted(((len(v), k) for k, v in by.items()), reverse=True)
    tot = len(tax); sp = len(by)
    print(f"{tot:,} records, {sp:,} species with a name")
    for cut in (1, 2, 5, 10, 50):
        k = sum(1 for c, _ in n if c >= cut)
        print(f"  species with >= {cut:>2} records: {k:>5,}  holding {sum(c for c, _ in n if c >= cut):>6,} records")
    print("  top:", ", ".join(f"{s} {c}" for c, s in n[:10]))
    jobs = [(i, by[s]) for c, s in n[:10] for i in by[s]]
    with Pool(8) as p:
        res = dict(p.map(nearest_same, jobs, chunksize=4))
    print("\nwithin-species nearest neighbour (edits over the shorter span), top 10 species")
    print(f"  {'species':30s} {'records':>7s} {'0':>5s} {'1':>5s} {'2':>5s} {'3-5':>5s} {'>5':>5s}  partial")
    for c, s in n[:10]:
        d = collections.Counter()
        for i in by[s]:
            x = res[i]; d["0" if x == 0 else "1" if x == 1 else "2" if x == 2 else "3-5" if x <= 5 else ">5"] += 1
        part = sum(1 for i in by[s] if tax[i]["grade"] not in ("complete", "near_complete"))
        print(f"  {s:30s} {c:>7} " + " ".join(f"{d[k]:>5}" for k in ("0", "1", "2", "3-5", ">5")) + f"  {part:>5}")
