"""cmp_top.py <acc dir> <TOP>: the reference classifier (candidate cap TOP) on amphid's own
sequences, vs amphid's calls; read-weighted."""
import csv, glob, os, subprocess, sys, collections
S = os.path.dirname(os.path.abspath(__file__))
d, top = sys.argv[1], sys.argv[2]
name = os.path.basename(d.rstrip("/"))
pre = f"{S}/amph.{name}"
rows = [r for p in glob.glob(d + "/*.clusters.tsv") for r in csv.DictReader(open(p), delimiter="\t")
        if r["origin"] != "chimera"]
uniq = {}
for r in rows:
    uniq.setdefault(r["consensus"], f"A{len(uniq)}")
with open(pre + ".asvs.tsv", "w") as fh:
    fh.write("asv\tsequence\n" + "".join(f"{a}\t{s}\n" for s, a in uniq.items()))
with open(pre + ".counts.tsv", "w") as fh:
    fh.write("sample\tasv\treads\n" + "".join(f"{r['sample_id']}\t{uniq[r['consensus']]}\t{r['reads']}\n" for r in rows))
out = f"{S}/reftop{top}.amph.{name}.tsv"
if os.path.exists(out): os.remove(out)
print(subprocess.run([sys.executable, f"{S}/ref_top.py", top, pre], capture_output=True, text=True).stdout.strip())
ref = {r["asv"]: r for r in csv.DictReader(open(out), delimiter="\t")}
tot = same = sp = 0; dis = collections.Counter()
for r in rows:
    n = int(r["reads"]); tot += n; x = ref[uniq[r["consensus"]]]
    a = r["rank"] if r["rank"] in ("species", "genus") else "other"
    b = x["gap_rank"] if x["gap_rank"] in ("species", "genus") else "other"
    sp += n * (a == "species")
    if a == b and (a == "other" or r["taxon"] == x["gap_name"]): same += n
    else: dis[(a, r["taxon"], r["gap"], r["next"], "| ref", b, x["gap_name"], x["gap"], x["next"])] += n
print(f"amphid species {100*sp/tot:.1f}%; same call as reference TOP={top} on {100*same/tot:.2f}% of reads")
for k, n in dis.most_common(6): print(" ", n, k)
