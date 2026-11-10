import csv, collections, sys
from multiprocessing import Pool
sys.path.insert(0, "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/bin")
from build import span_distance, read_fasta
B = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/work/build_v0.1.0"
seqs = read_fasta(B + "/refdb.fasta")
tax = {r["id"]: r for r in csv.DictReader(open(B + "/taxonomy.tsv"), delimiter="\t")}
mem = collections.defaultdict(list)
for r in csv.DictReader(open(B + "/members.tsv"), delimiter="\t"):
    if r["role"] == "kept": mem[r["id"]].append(r["accession"])
hc = [i for i, t in tax.items() if t["species"] == "Haemonchus contortus"]
hp = [i for i, t in tax.items() if t["species"] == "Haemonchus placei"]
def near(i):
    return i, min((span_distance(seqs[i], seqs[j]), j) for j in hp)
if __name__ == "__main__":
    with Pool(8) as p: res = dict(p.map(near, hc, chunksize=4))
    d = collections.Counter(min(v[0], 9) for v in res.values())
    print("H. contortus records by edits to the nearest H. placei record (9 = 9+):", dict(sorted(d.items())))
    close = [i for i in hc if res[i][0] <= 3]
    print(f"{len(close)} H. contortus records within 3 edits of H. placei; grades:",
          dict(collections.Counter(tax[i]["grade"] for i in close)),
          "; accessions per record:", dict(collections.Counter(len(mem[i]) for i in close)))
    for i in sorted(close, key=lambda i: res[i][0])[:10]:
        print(f"  {i} {res[i][0]} edits  {tax[i]['grade']:13s} len {len(seqs[i])}  {','.join(mem[i][:3])}{' +' + str(len(mem[i])-3) if len(mem[i])>3 else ''}")
    for i in ("ITS2_299134f6aa91", "ITS2_d2974ff6ffa8"):
        print("  decisive record", i, tax[i]["grade"], len(seqs[i]), "edits to placei", res[i][0], mem[i][:4])
    top = lambda ids: max(ids, key=lambda i: len(mem[i]))
    a, b = top(hc), top(hp)
    print("most-supported H. contortus vs H. placei record:", span_distance(seqs[a], seqs[b]), "edits",
          f"({len(mem[a])} and {len(mem[b])} accessions)")
