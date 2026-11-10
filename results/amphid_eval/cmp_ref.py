"""Run nemabiome-its2-refdb's classify_asvs.py on amphid's sequences and
compare its gap-rule calls with amphid's, read-weighted.
usage: cmp_ref.py <clusters.tsv>... (writes <scratch>/cmp.asvs.tsv)"""
import csv, subprocess, sys, os, collections
R = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb"
S = os.path.dirname(os.path.abspath(__file__))
rows = []
for p in sys.argv[1:]:
    rows += [r for r in csv.DictReader(open(p), delimiter="\t") if r["origin"] != "chimera"]
uniq = {}
for r in rows:
    uniq.setdefault(r["consensus"], f"A{len(uniq)}")
with open(f"{S}/cmp.asvs.tsv", "w") as fh:
    fh.write("asv\tsequence\n")
    for s, a in uniq.items():
        fh.write(f"{a}\t{s}\n")
subprocess.run([sys.executable, f"{R}/bin/classify_asvs.py", f"{S}/cmp.asvs.tsv",
                f"{R}/work/build_v0.1.0", f"{S}/cmp.out.tsv", "--workers", "8"], check=True)
ref = {r["asv"]: r for r in csv.DictReader(open(f"{S}/cmp.out.tsv"), delimiter="\t")}
tot = agree = 0
dis = collections.Counter()
sp_a = sp_r = 0
for r in rows:
    n = int(r["reads"]); tot += n
    x = ref[uniq[r["consensus"]]]
    a = (r["rank"] if r["rank"] in ("species", "genus", "family") else "other", r["taxon"])
    b = (x["gap_rank"] if x["gap_rank"] in ("species", "genus", "family") else "other", x["gap_name"])
    sp_a += n * (a[0] == "species"); sp_r += n * (b[0] == "species")
    if a == b or (a[0] == b[0] == "other"):
        agree += n
    else:
        dis[(a, b, r["identity"], x["identity"], r["gap"], x["gap"], r["next"], x["next"])] += n
print(f"reads {tot}: same call {100*agree/tot:.2f}%; species amphid {100*sp_a/tot:.2f}% ref {100*sp_r/tot:.2f}%")
for k, n in dis.most_common(15):
    print(n, k)
