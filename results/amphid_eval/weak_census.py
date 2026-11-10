"""Weak records (partial or unmeasured, one accession) by whether their species also has robust ones,
and how far each weak record sits from its species' nearest robust record."""
import collections, csv, sys
from multiprocessing import Pool
sys.path.insert(0, "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/bin")
from build import span_distance, read_fasta
B = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/work/build_v0.1.0"
seqs = read_fasta(B + "/refdb.fasta")
tax = {r["id"]: r for r in csv.DictReader(open(B + "/taxonomy.tsv"), delimiter="\t")}
nm = collections.Counter(r["id"] for r in csv.DictReader(open(B + "/members.tsv"), delimiter="\t") if r["role"] == "kept")
FULL = ("complete", "near_complete")
weak = {i for i, t in tax.items() if t["grade"] not in FULL and nm[i] == 1}
full = {i for i, t in tax.items() if t["grade"] in FULL}
multi = {i for i in tax if nm[i] >= 2}
by = collections.defaultdict(list)
for i, t in tax.items():
    by[t["species"]].append(i)
def nearest(i):
    rob = [j for j in by[tax[i]["species"]] if j not in weak]
    return i, min(span_distance(seqs[i], seqs[j]) for j in rob)
if __name__ == "__main__":
    named = [i for i in weak if tax[i]["species"]]
    print(f"{len(weak):,} weak records: {len(named):,} name a species, {len(weak) - len(named):,} do not (shared/unnamed/genus-level)")
    a = [i for i in named if any(j not in weak for j in by[tax[i]["species"]])]
    b = [i for i in named if any(j in full and j in multi for j in by[tax[i]["species"]])]
    only = [i for i in named if i not in a]
    print(f"  in species that also have a robust record (full-length or >= 2 accessions): {len(a):,} "
          f"({len({tax[i]['species'] for i in a}):,} species)")
    print(f"  in species with a full-length, multi-accession record:                      {len(b):,} "
          f"({len({tax[i]['species'] for i in b}):,} species)")
    print(f"  the only evidence for their species:                                        {len(only):,} "
          f"({len({tax[i]['species'] for i in only}):,} species)")
    top = collections.Counter(tax[i]["species"] for i in a).most_common(10)
    print("  most weak records, among species with robust ones:", ", ".join(f"{s} {n}" for s, n in top))
    with Pool(8) as p:
        d = dict(p.map(nearest, a, chunksize=20))
    c = collections.Counter("0" if v == 0 else "1" if v == 1 else "2" if v == 2 else "3-5" if v <= 5 else ">5" for v in d.values())
    print("  edits from each to its species' nearest robust record:", {k: c[k] for k in ("0", "1", "2", "3-5", ">5")})
    hc = [i for i in a if tax[i]["species"] == "Haemonchus contortus"]
    ch = collections.Counter("0" if d[i] == 0 else "1" if d[i] == 1 else "2" if d[i] == 2 else "3-5" if d[i] <= 5 else ">5" for i in hc)
    print(f"  H. contortus: {len(hc)} weak of {len(by['Haemonchus contortus'])}; edits to nearest robust:",
          {k: ch[k] for k in ("0", "1", "2", "3-5", ">5")})
