#!/usr/bin/env python3
"""HTML report of one project's real-data runs (bin/real_data.sh).

    bin/real_data_report.py <out dir> <project> <title> <report.html> [notes.html]

Reads <out dir>/<project>/{250,150}/trim_stats.tsv, every DADA2 run
<out dir>/<project>.<set> (track, counts, and each arm's calls, as
bin/compare_arms.py does) and work/testdata/<project>.runs.tsv for the
host of each run (from work/testdata/<project>.hosts.tsv, run_accession
and host, where the run titles do not name it). Where
work/testdata/<project>.replicates.tsv (run_accession, sample) pairs
libraries of one sample, the page can also count by sample, keeping only
ASVs found in every library of it. An ASV no record aligns to (blank
`nearest` in the gap rule's table; off-target products, primer dimers)
is shown as "not ITS2-like", apart from unassigned ("unplaced ITS2"), when a
classifier leaves it above genus. One self-contained page: trimming, the 2x250 merge
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
        ("ours.nemits2_dev", "tiered", "Tiered rule (gap + tentative), this database"),
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
    """A rule's calls, or None if the table predates it; a tentative
    species (tiered rule) is labelled so."""
    rows = tsv(f)
    if rows and rule + "_rank" not in rows[0]:
        return None
    return {r["asv"]: (r[rule + "_rank"], r[rule + "_name"]
                       + (" (tentative)" if rule == "tiered" and r.get("tentative") == "yes" else ""))
            for r in rows}


NONE = "not ITS2-like"


def call(rank, name, stray=False):
    """(level, genus, label): level species / genus / higher (family or
    above) / unassigned (unplaced ITS2), or none for a stray ASV (not ITS2-like) left
    above genus."""
    if stray and rank not in ("species", "genus"):
        return "none", "", NONE
    if rank == "species":
        return "species", name.split()[0], name
    if rank == "genus":
        return "genus", name, name + " (genus)"
    if rank == "unassigned" or not name:
        return "unassigned", "", "unplaced ITS2"
    return "higher", "", f"{name} ({rank})"


runs = {r["run_accession"]: r for r in tsv(f"work/testdata/{P}.runs.tsv")}
if os.path.exists(f"work/testdata/{P}.hosts.tsv"):     # run_accession, host
    host = {r["run_accession"]: r["host"] for r in tsv(f"work/testdata/{P}.hosts.tsv")}
else:
    host = {k: re.search(r"(cattle|sheep|goat|pig|kangaroo|rabbit)", r["experiment_title"]).group(1)
            for k, r in runs.items()}
alias = {k: r["library_name"] or r["sample_alias"] for k, r in runs.items()}
trim_dir = f"{OUT}/{P}/250" if os.path.exists(f"{OUT}/{P}/250/trim_stats.tsv") else f"{OUT}/{P}/150"
trim = {r["run"]: r for r in tsv(f"{trim_dir}/trim_stats.tsv")}
rep_f = f"work/testdata/{P}.replicates.tsv"
group = {r["run_accession"]: r["sample"] for r in tsv(rep_f)} if os.path.exists(rep_f) else {}

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
             arms={}, by_sample={})
    stray = {r["asv"] for r in tsv(f"{pre}.ours.nemits2_dev.tsv")
             if r["gap_rank"] == "unassigned" and not r["nearest"]} \
        if os.path.exists(f"{pre}.ours.nemits2_dev.tsv") else set()
    pooled = collections.defaultdict(collections.Counter)    # sample -> ASVs in all its libraries
    libs = collections.defaultdict(list)
    for smp in counts:
        if smp in group:
            libs[group[smp]].append(smp)
    for g, ls in libs.items():
        keep = set.intersection(*({a for a, n in counts[l].items() if n} for l in ls))
        for l in ls:
            pooled[g].update({a: n for a, n in counts[l].items() if a in keep})
    for suf, rule, alab in ARMS:
        f = f"{pre}.{suf}.tsv"
        if not os.path.exists(f):
            continue
        calls = ours(f, rule) if rule else kmer(f)
        if calls is None:
            continue
        def comp_of(c):
            comp = collections.Counter()
            for asv, n in c.items():
                comp[call(*calls.get(asv, ("unassigned", "")), asv in stray)] += n
            return [[lv, g, lab, n] for (lv, g, lab), n in comp.most_common()]
        s["arms"][alab] = {smp: comp_of(c) for smp, c in counts.items()}
        if pooled:
            s["by_sample"][alab] = {g: comp_of(c) for g, c in pooled.items()}
    data["sets"].append(s)

samples = sorted(trim, key=lambda k: (host[k], alias[k]))
for k in samples:
    t = trim[k]
    data["samples"].append(dict(run=k, alias=alias[k], host=host[k], pairs=int(t["pairs"]),
                                kept=int(t["kept"]), nc2=int(t["read_2_at_NC2"]),
                                in28s=int(t["read_2_in_28S"])))

if group:
    data["groups"] = [dict(id=g, host=host[ls[0]], runs=ls) for g, ls in sorted(
        {g: sorted(r for r in group if group[r] == g and r in trim) for g in set(group.values())}.items(),
        key=lambda x: (host[x[1][0]] if x[1] else "", x[0])) if ls]

page = open(os.path.join(os.path.dirname(__file__), "real_data_report.html")).read()
page = (page.replace("{{TITLE}}", html.escape(TITLE)).replace("{{NOTES}}", NOTES)
        .replace("{{DATA}}", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))
os.makedirs(os.path.dirname(DEST) or ".", exist_ok=True)
open(DEST, "w").write(page)
print(f"{DEST}: {len(data['sets'])} sets, {len(samples)} samples, {os.path.getsize(DEST):,} bytes")
