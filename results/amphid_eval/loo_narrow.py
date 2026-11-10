"""Leave-one-out (as bin/its2_leave_one_out.py; v0.1.0, no candidate cap), amphid's gap rule
as released vs with single-accession partial records not counted as another species."""
import collections, csv, os, random, subprocess, sys, tempfile
from multiprocessing import Pool
R = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb"
A = "/Users/sd21/Documents/WORK/informatics_workbook/nf-amphid-pipeline/bin"
B = R + "/work/build_v0.1.0"
sys.path[:0] = [R + "/bin", A]
G = {}
def init():
    import ampseqlib as al, assign_db as adb
    G["adb"] = adb
    G["seqs"] = dict(al.read_fasta(B + "/refdb.fasta"))
    G["tax"] = tax = al.load_taxonomy(B + "/taxonomy.tsv")
    nm = collections.Counter(r["id"] for r in csv.DictReader(open(B + "/members.tsv"), delimiter="\t") if r["role"] == "kept")
    import json
    G["sets"] = {k: set(v) for k, v in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "narrow_sets.json"))).items()}
def outcome(c, t, tax, adb):
    if c["rank"] == "species":
        return "right species" if tax[c["best_hit"]]["species"] == t["species"] else "WRONG species"
    if c["rank"] == "genus":
        return "genus, right" if tax[c["best_hit"]]["genus"] == t["genus"] else "WRONG genus"
    return "family or higher" if c["depth"] >= 0 and c["identity"] >= 0.95 else "unassigned"
def job(j):
    (qi, mode, view), x, y, cs = j
    adb, tax, seqs = G["adb"], G["tax"], G["seqs"]
    parts = [(x + "N" * 47, "whole")] if view == "whole" else [(x, "start"), (y + "N" * 47, "end")]
    sc = adb.span_scores(parts, cs, seqs, 47)
    a = adb.span_call(sc, tax, "gap", adb.GAP_FLOORS)
    best_sp = {tax[c]["species"] for c in a["ties"]}
    t = tax[qi]
    out = [outcome(a, t, tax, adb)]
    for k in ("suspect", "suspect+redundant"):
        drop = G["sets"][k]
        b = adb.span_call([s for s in sc if s[2] not in drop or tax[s[2]]["species"] in best_sp], tax, "gap", adb.GAP_FLOORS)
        out.append(outcome(b, t, tax, adb))
    return ((qi, mode, view), *out)
def main():
    init()
    seqs, tax = G["seqs"], G["tax"]
    holds = collections.defaultdict(set)
    for r in csv.DictReader(open(B + "/members.tsv"), delimiter="\t"):
        if r["role"] == "kept" and r["id"]: holds[r["label"].split(";")[-1]].add(r["id"])
    rt = {r["id"]: r for r in csv.DictReader(open(B + "/taxonomy.tsv"), delimiter="\t")}
    by = collections.defaultdict(list)
    for i, t in rt.items():
        if t["deepest_rank"] == "species" and t["grade"] in ("complete", "near_complete"): by[t["species"]].append(i)
    rng = random.Random(7)
    queries = [q for sp in sorted(by) for q in rng.sample(sorted(by[sp]), min(5, len(by[sp])))]
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "q.fa"); open(f, "w").writelines(f">{i}\n{seqs[i]}\n" for i in queries)
        out = subprocess.run(["minimap2", "-t", "8", "-k", "11", "-w", "5", "-f", "1000000", "-N", "5000", "-p", "0.3",
                              B + "/refdb.fasta", f], capture_output=True, text=True, check=True).stdout
    cand = collections.defaultdict(dict)
    for line in out.splitlines():
        x = line.split("\t"); cand[x[0]][x[5]] = max(cand[x[0]].get(x[5], 0), int(x[9]))
    jobs = []
    for qi in queries:
        sp = rt[qi]["species"]; ranked = sorted(cand[qi], key=lambda c: -cand[qi][c])
        for mode, drop in (("hidden", holds[sp] | {qi}), ("present", {qi})):
            if mode == "present" and not (holds[sp] - {qi}): continue
            cs = [c for c in ranked if c not in drop]; s = seqs[qi]
            jobs.append(((qi, mode, "whole"), s, "", cs))
            if len(s) > 213: jobs.append(((qi, mode, "2x150"), s[:130], s[-83:], cs))
    with Pool(8, initializer=init) as p:
        res = p.map(job, jobs, chunksize=50)
    cats = ["right species", "WRONG species", "genus, right", "WRONG genus", "family or higher", "unassigned"]
    print(f"{len(queries):,} queries, {len(by):,} species; " + ", ".join(f"{k}: {len(v):,} records" for k, v in G["sets"].items()))
    for mode in ("hidden", "present"):
        for view in ("whole", "2x150"):
            sel = [r for r in res if r[0][1] == mode and r[0][2] == view]
            print(f"== {mode}, {view} ({len(sel):,})\n   {'':22s}" + "".join(f"{c:>17s}" for c in cats))
            for name, i in (("as released", 1), ("suspect", 2), ("suspect+redundant", 3)):
                c = collections.Counter(r[i] for r in sel)
                print(f"   {name:22s}" + "".join(f"{100*c[k]/len(sel):16.1f}%" for k in cats))
            # which wrong calls are new
            for name, i in (("suspect", 2), ("suspect+redundant", 3)):
                new = collections.Counter(tax[r[0][0]]["genus"] for r in sel if r[i] == "WRONG species" and r[1] != "WRONG species")
                if new: print(f"   new wrong-species calls, {name}:", ", ".join(f"{g} {n}" for g, n in new.most_common(8)))
if __name__ == "__main__":
    main()
