"""V4 regression: main's assign_db(_pairs).py vs the branch's on the fixture, outputs diffed."""
import gzip, json, os, random, subprocess, sys, filecmp, glob
A = "/Users/sd21/Documents/WORK/informatics_workbook/nf-amphid-pipeline"
S = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [f"{A}/tests", f"{A}/bin"]
import make_db_fixture as fx, simulate as sim, ampseqlib as al
W = f"{S}/v4"; os.makedirs(W, exist_ok=True)
rng = random.Random(20261028)
species = fx.build_tree(rng)
rel = f"{W}/release"
if not os.path.exists(rel):
    fx.write_release(rel, species)
fx.write_primers(f"{W}/primers.tsv")
model = sim.ERROR_MODELS["nanopore_hac"]
sheet = [{"sample_id": s} for s, _t, _c in fx.SAMPLES]
cfg = dict(fx.sim.DEFAULTS); cfg.update(al.load_assay(fx.ASSAY))
fx.write_illumina(f"{W}/run", sheet, species, cfg, random.Random(9))
for sid, _t, comp in fx.SAMPLES:
    with gzip.open(f"{W}/{sid}.fastq.gz", "wt") as fh:
        k = 0
        for src, n in comp.items():
            for _ in range(n):
                m = fx.amplicon(species[src]["seq"])
                if rng.random() < 0.5: m = al.revcomp(m)
                m, _s, _i = sim.apply_errors(rng, m, model)
                fh.write(f"@{src}|{k}\n{m}\n+\n{'5' * len(m)}\n"); k += 1
for tag, b in (("main", f"{S}/orig_bin"), ("branch", f"{A}/bin")):
    for sid, _t, _c in fx.SAMPLES:
        subprocess.run([sys.executable, f"{b}/assign_db.py", "--reads", f"{W}/{sid}.fastq.gz", "--database", rel,
                        "--primers", f"{W}/primers.tsv", "--sample-id", sid, "--outdir", f"{W}/{tag}/ont"],
                       check=True, capture_output=True)
        for rl in ("150", "300"):
            subprocess.run([sys.executable, f"{b}/assign_db_pairs.py", "--r1", f"{W}/run/{sid}_R1.fastq.gz",
                            "--r2", f"{W}/run/{sid}_R2.fastq.gz", "--database", rel, "--primers", f"{W}/primers.tsv",
                            "--assay", fx.ASSAY, "--sample-id", sid, "--read-length", rl,
                            "--outdir", f"{W}/{tag}/il{rl}"], check=True, capture_output=True)
same = diff = 0
for f in sorted(glob.glob(f"{W}/main/*/*")):
    g = f.replace("/main/", "/branch/")
    if f.endswith(".json"):
        a, b = json.load(open(f)), json.load(open(g))
        for x in (a, b): x.pop("elapsed_s", None); x["database"].pop("path", None)
        ok = a == b
    else:
        ok = filecmp.cmp(f, g, shallow=False)
    same += ok; diff += not ok
    if not ok: print("DIFFERS:", f)
print(f"{same} files identical, {diff} differ")
