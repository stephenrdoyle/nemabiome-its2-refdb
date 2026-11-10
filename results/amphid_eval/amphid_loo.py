"""amphid_loo.py <build> <out prefix> [--cands ref|amphid]: nemabiome-its2-refdb's leave-one-out
(bin/its2_leave_one_out.py: same queries, same hidden/present removal, same views) scored by
amphid's span scoring and gap rule (assign_db.span_scores / span_call), side by side with the
reference's _score/call on the same jobs."""
import argparse, collections, csv, os, random, subprocess, sys, tempfile
from multiprocessing import Pool
R = "/Users/sd21/Documents/WORK/informatics_workbook/nemabiome-its2-refdb/bin"
A = "/Users/sd21/Documents/WORK/informatics_workbook/nf-amphid-pipeline/bin"
sys.path[:0] = [R, A]
G = {}

def init(build):
    import its2_classify as ic, ampseqlib as al, assign_db as adb
    from build import read_fasta
    seqs = read_fasta(os.path.join(build, "refdb.fasta"))
    G.update(seqs=seqs, rtax=ic.load_taxonomy(build), atax=al.load_taxonomy(os.path.join(build, "taxonomy.tsv")),
             ic=ic, adb=adb)
    ic._init(seqs, G["rtax"])

def job(j):
    key, view, x, y, cs = j
    ic, adb = G["ic"], G["adb"]
    _k, d, sp, ties, gap, cv, _n = ic._score((key, "whole" if view == "whole" else "window",
                                              x + ("N" * 47 if view == "whole" else ""),
                                              (y + "N" * 47) if y else "", cs))
    rc = ic.call(G["rtax"], ties, d, sp, "gap", gap, cv)
    parts = [(x + "N" * 47, "whole")] if view == "whole" else [(x, "start"), (y + "N" * 47, "end")]
    c = adb.span_call(adb.span_scores(parts, cs, G["seqs"], 47), G["atax"], "gap", adb.GAP_FLOORS)
    lin = G["atax"].get(c["best_hit"], {})
    # Below the family floor the reference makes no call; amphid keeps a
    # rank but labels it "unclassified, NN% to X". Both are unassigned.
    ac = (c["rank"] if c["rank"] in ("species", "genus") else
          ("family or higher" if c["depth"] >= 0 and c["identity"] >= adb.GAP_FLOORS["family"] else ""),
          lin.get(c["rank"], "") if c["rank"] else "")
    return key, rc, ac

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("build"); ap.add_argument("prefix")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cap", type=int, default=60, help="candidates per query; 0 = all")
    ap.add_argument("-N", type=int, default=300, help="minimap2 -N")
    a = ap.parse_args()
    from build import read_fasta
    seqs = read_fasta(os.path.join(a.build, "refdb.fasta"))
    tax = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.build, "taxonomy.tsv")), delimiter="\t")}
    holds = collections.defaultdict(set)
    for r in csv.DictReader(open(os.path.join(a.build, "members.tsv")), delimiter="\t"):
        if r["role"] == "kept" and r["id"]:
            holds[r["label"].split(";")[-1]].add(r["id"])
    by_sp = collections.defaultdict(list)
    for i, t in tax.items():
        if t["deepest_rank"] == "species" and t["grade"] in ("complete", "near_complete"):
            by_sp[t["species"]].append(i)
    rng = random.Random(7)
    queries = [q for sp in sorted(by_sp) for q in rng.sample(sorted(by_sp[sp]), min(5, len(by_sp[sp])))]
    with tempfile.TemporaryDirectory() as d:
        q = os.path.join(d, "q.fa")
        with open(q, "w") as fh:
            fh.writelines(f">{i}\n{seqs[i]}\n" for i in queries)
        out = subprocess.run(["minimap2", "-t", str(a.workers), "-k", "11", "-w", "5", "-f", "1000000",
                              "-N", str(a.N), "-p", "0.3", os.path.join(a.build, "refdb.fasta"), q],
                             capture_output=True, text=True, check=True).stdout
    cand = collections.defaultdict(dict)
    for line in out.splitlines():
        f = line.split("\t")
        cand[f[0]][f[5]] = max(cand[f[0]].get(f[5], 0), int(f[9]))
    jobs = []
    for qi in queries:
        sp = tax[qi]["species"]
        ranked = sorted(cand[qi], key=lambda c: -cand[qi][c])
        for mode, drop in (("hidden", holds[sp] | {qi}), ("present", {qi})):
            if mode == "present" and not (holds[sp] - {qi}):
                continue
            cs = [c for c in ranked if c not in drop]
            if a.cap: cs = cs[:a.cap]
            s = seqs[qi]
            jobs.append(((qi, mode, "whole"), "whole", s, "", cs))
            if len(s) > 130 + 83:
                jobs.append(((qi, mode, "2x150"), "2x150", s[:130], s[-83:], cs))
    with Pool(a.workers, initializer=init, initargs=(a.build,)) as pool:
        res = list(pool.imap_unordered(job, jobs, chunksize=50))
    def outcome(qi, c):
        rank, name = c; t = tax[qi]
        if rank == "species": return "right species" if name == t["species"] else "WRONG species"
        if rank == "genus": return "genus, right" if name == t["genus"] else "WRONG genus"
        return "family or higher" if rank else "unassigned"
    cats = ["right species", "WRONG species", "genus, right", "WRONG genus", "family or higher", "unassigned"]
    lines = [f"{len(queries):,} queries from {len(by_sp):,} species, {a.build}; candidates: minimap2 -N {a.N}, cap {a.cap or 'none'}"]
    nc = [len(j[4]) for j in jobs]; lines.append(f"candidates per job: median {sorted(nc)[len(nc)//2]}, max {max(nc)}")
    same = sum(outcome(k[0], rc) == outcome(k[0], ac) for k, rc, ac in res)
    lines.append(f"amphid and reference give the same outcome on {same:,} of {len(res):,} jobs")
    for mode in ("hidden", "present"):
        for view in ("whole", "2x150"):
            sel = [x for x in res if x[0][1] == mode and x[0][2] == view]
            lines.append(f"== {mode}, {view} ({len(sel):,})\n   {'':10s}" + "".join(f"{c:>17s}" for c in cats))
            for name, i in (("reference", 1), ("amphid", 2)):
                c = collections.Counter(outcome(x[0][0], x[i]) for x in sel)
                lines.append(f"   {name:10s}" + "".join(f"{100 * c[k] / len(sel):16.1f}%" for k in cats))
    open(a.prefix + ".txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

if __name__ == "__main__":
    main()
