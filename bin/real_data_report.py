#!/usr/bin/env python3
"""HTML report of one project's real-data runs (bin/real_data.sh).

    bin/real_data_report.py <out dir> <project> <title> <report.html> [notes.html]

Reads <out dir>/<project>/{250,150}/trim_stats.tsv, every DADA2 run
<out dir>/<project>.<set> (track, counts, and each arm's calls, as
bin/compare_arms.py does) and work/testdata/<project>.runs.tsv for the
host of each run. One self-contained page: trimming, the 2x250 merge
options against 2x150, the share of reads each classifier names to
species, genus or above, per-sample composition, and a per-sample table.
[notes.html] is a fragment placed under the title (the findings).
"""
import csv, collections, glob, html, json, os, re, sys

OUT, P, TITLE, DEST = sys.argv[1:5]
NOTES = open(sys.argv[5]).read() if len(sys.argv) > 5 else ""
RANKS = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
SETS = [("150", "2x150, concatenated"), ("250.merge", "2x250, merged"),
        ("250.concat", "2x250, concatenated"), ("250.hybrid", "2x250, hybrid")]
ARMS = [("ours.nemits2_dev", "gap", "Gap rule, this database"),
        ("idtaxa.nemits2_dev.trim28s", None, "IdTaxa, this database, 28S cut"),
        ("tax.nemits2_dev.trim28s", None, "assignTaxonomy, this database, 28S cut"),
        ("idtaxa.nemabiome_20260317", None, "IdTaxa, nemabiome backbone"),
        ("tax.nemabiome_20260317", None, "assignTaxonomy, nemabiome backbone")]


def tsv(f):
    return list(csv.DictReader(open(f), delimiter="\t"))


def kmer(f):
    """assignTaxonomy / IdTaxa: deepest named rank, columns by position."""
    out = {}
    for r in tsv(f):
        cols = [c for c in r if not c.startswith("boot") and c != "asv"][:8]
        rank, name = "unassigned", ""
        for i, c in enumerate(cols):
            v = r[c]
            if v and v != "NA":
                rank, name = RANKS[i], v.replace("_", " ")
            else:
                break
        out[r["asv"]] = (rank, name)
    return out


def ours(f, rule):
    return {r["asv"]: (r[rule + "_rank"], r[rule + "_name"]) for r in tsv(f)}


def call(rank, name):
    """(level, genus, label): level species / genus / higher."""
    if rank == "species":
        return "species", name.split()[0], name
    if rank == "genus":
        return "genus", name, name + " (genus)"
    return "higher", "", (f"{name} ({rank})" if name else "unassigned")


runs = {r["run_accession"]: r for r in tsv(f"work/testdata/{P}.runs.tsv")}
host = {k: re.search(r"(cattle|sheep|goat|pig|kangaroo|rabbit)", r["experiment_title"]).group(1)
        for k, r in runs.items()}
alias = {k: r["library_name"] or r["sample_alias"] for k, r in runs.items()}
trim = {r["run"]: r for r in tsv(f"{OUT}/{P}/250/trim_stats.tsv")}

data = dict(project=P, title=TITLE, sets=[], arms=[a[2] for a in ARMS], samples=[])
for sid, slab in SETS:
    pre = f"{OUT}/{P}.{sid}"
    if not os.path.exists(pre + ".counts.tsv"):
        continue
    counts = collections.defaultdict(collections.Counter)
    for r in tsv(pre + ".counts.tsv"):
        counts[r["sample"]][r["asv"]] += int(r["reads"])
    track = {r["sample"]: r for r in tsv(pre + ".track.tsv")}
    s = dict(id=sid, label=slab, asvs=len({a for c in counts.values() for a in c}),
             track={k: {c: int(v[c]) for c in ("input", "used", "nonchim")} for k, v in track.items()},
             arms={})
    for suf, rule, alab in ARMS:
        f = f"{pre}.{suf}.tsv"
        if not os.path.exists(f):
            continue
        calls = ours(f, rule) if rule else kmer(f)
        per = {}
        for smp, c in counts.items():
            comp = collections.Counter()
            for asv, n in c.items():
                comp[call(*calls.get(asv, ("unassigned", "")))] += n
            per[smp] = [[lv, g, lab, n] for (lv, g, lab), n in comp.most_common()]
        s["arms"][alab] = per
    data["sets"].append(s)

samples = sorted(trim, key=lambda k: (host[k], alias[k]))
for k in samples:
    t = trim[k]
    data["samples"].append(dict(run=k, alias=alias[k], host=host[k], pairs=int(t["pairs"]),
                                kept=int(t["kept"]), nc2=int(t["read_2_at_NC2"]),
                                in28s=int(t["read_2_in_28S"])))

page = open(os.path.join(os.path.dirname(__file__), "real_data_report.html")).read()
page = (page.replace("{{TITLE}}", html.escape(TITLE)).replace("{{NOTES}}", NOTES)
        .replace("{{DATA}}", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))
os.makedirs(os.path.dirname(DEST) or ".", exist_ok=True)
open(DEST, "w").write(page)
print(f"{DEST}: {len(data['sets'])} sets, {len(samples)} samples, {os.path.getsize(DEST):,} bytes")
