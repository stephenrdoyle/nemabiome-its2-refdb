#!/usr/bin/env python3
"""Score classifier arms on the simulated nemabiome mocks against their truth.

    bin/score_mocks.py <mock dir> <source fasta> <overrides.tsv> <backbone seqs fasta>

The mocks (Jesudoss Chelladurai et al. 2025, Parasites & Vectors;
Zenodo 10.5281/zenodo.17019940) are ART-Illumina 2x250 reads simulated
from nemabiome v1.0.0 records; each read is named for its species, not
its record. <mock dir> holds one bin/dada2_nemabiome.R run (<name>.asvs.tsv,
.counts.tsv) and the arms' calls, as for bin/compare_arms.py.

An ASV's truth is the species of the source records it equals (or sits
inside, or holds whole). An ASV matching none is mostly a source copy
carrying simulated sequencing errors DADA2 kept as a variant: if it is
<= NEAR edits from sources of one species only, it is "noisy" with that
species as truth, and scored in a table of its own - the closest the mocks
come to a read whose exact variant is not in the database. Otherwise
"no source". Two kinds of source record
are scored apart too, because their mock label is not the species:
  curated   records this database drops as mislabelled
            (curation/overrides.tsv, found by exact sequence in the backbone)
  synonym   Cylicostephanus poculatus = Petrovinema poculatum: the mock
            lists both; one species
Per arm, read-weighted, at species and at genus: recall (right calls /
reads) and precision (right calls / reads called at that rank or deeper;
a species call is right at genus when its genus is right). Every mock
species is in the database, so there is no over-classification to
measure here (bin/loo_metrics.py does, by leaving species out). And per
sample, the species
abundance error: half the summed absolute difference between called and
true species shares (reads not named to species count as error) - over
clean ASVs ("abund. err", the classifier alone), and over all reads
against each sample's composition from the read names ("vs mock", which
also carries what DADA2 lost, the same for every arm).
"""

from __future__ import annotations

import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import RANKS, infix_distance, read_fasta  # noqa: E402

SYNONYMS = {"Cylicostephanus poculatus": "Petrovinema poculatum",
            "Parapoteriostomum mettami": "Cylicodontophorus mettami"}
NEAR = 8        # an ASV this many edits or fewer from a source is a noisy copy of it


def canon(sp):
    return SYNONYMS.get(sp, sp)


def k_mer_calls(path):
    """assignTaxonomy / IdTaxa table: (rank, name) by column position."""
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        cols = [c for c in r if not c.startswith("boot") and c != "asv"][:8]
        rank, name = "", ""
        for i, c in enumerate(cols):
            if r[c] and r[c] != "NA":
                rank, name = RANKS[i], r[c].replace("_", " ")
            else:
                break
        out[r["asv"]] = (rank, name)
    return out


def our_calls(path):
    return {r["asv"]: (r["gap_rank"] if r["gap_rank"] != "unassigned" else "", r["gap_name"])
            for r in csv.DictReader(open(path), delimiter="\t")}


def main() -> int:
    d, src_fa, overrides, backbone = sys.argv[1:5]
    name = os.path.basename(d.rstrip("/"))
    P = os.path.join(d, name)
    # sources: sequence -> species; curated sources -> what they really are
    src = []
    for line in open(src_fa):
        line = line.strip()
        if line.startswith(">"):
            src.append([canon(line[1:].replace("_", " ")), ""])
        elif line:
            src[-1][1] += line.upper()
    bb = collections.defaultdict(list)
    for acc, s in read_fasta(backbone).items():
        bb[s.upper()].append(acc.split(".")[0])
    dropped = {r[0].split(".")[0]: r[3] for r in csv.reader(open(overrides), delimiter="\t")
               if r and not r[0].startswith("#") and r[1] == "drop"}
    curated = {s for _, s in src if any(a in dropped for a in bb.get(s, []))}

    asvs = {r["asv"]: r["sequence"].upper() for r in csv.DictReader(open(P + ".asvs.tsv"), delimiter="\t")}
    reads = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(P + ".counts.tsv"), delimiter="\t"):
        reads[r["sample"]][r["asv"]] += int(r["reads"])

    truth, kind = {}, {}
    for a, s in asvs.items():
        hit = [(sp, t) for sp, t in src if s == t] or [(sp, t) for sp, t in src if s in t or t in s]
        sps = {sp for sp, _ in hit}
        if not hit:
            dist = [(min(infix_distance(t, s), infix_distance(s, t)), sp) for sp, t in src]
            best = min(x for x, _ in dist)
            near = {sp for x, sp in dist if x == best}
            kind[a] = "noisy" if best <= NEAR and len(near) == 1 else "no source"
            truth[a] = near.pop() if kind[a] == "noisy" else ""
            continue
        elif any(t in curated for _, t in hit):
            kind[a] = "curated"
        elif len(sps) > 1:
            kind[a] = "ambiguous source"
        else:
            kind[a] = "clean"
        truth[a] = sps.pop() if len(sps) == 1 else ""

    arms = {}
    for label, f, reader in (
            ("assignTaxonomy x nemabiome", ".tax.nemabiome_20260317.tsv", k_mer_calls),
            ("IdTaxa x nemabiome", ".idtaxa.nemabiome_20260317.tsv", k_mer_calls),
            ("assignTaxonomy x this db", ".tax.nemits2_dev.trim28s.tsv", k_mer_calls),
            ("IdTaxa x this db", ".idtaxa.nemits2_dev.trim28s.tsv", k_mer_calls),
            ("gap rule x this db", ".ours.nemits2_dev.tsv", our_calls)):
        if os.path.exists(P + f):
            arms[label] = collections.defaultdict(lambda: ("", ""), reader(P + f))

    def outcome(call, t):
        rank, nm = call
        nm = canon(nm)
        if rank == "species":
            return "right species" if nm == t else "WRONG species"
        if rank == "genus":
            return "genus, right" if nm == t.split()[0] else "WRONG genus"
        return "higher/unassigned"
    cats = ["right species", "WRONG species", "genus, right", "WRONG genus", "higher/unassigned"]
    metrics = ["species recall", "precision", "genus recall", "precision"]

    def call_genus(call):
        rank, nm = call
        return canon(nm).split()[0] if rank == "species" else nm

    mock = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(os.path.join(d, "truth.tsv")), delimiter="\t"):
        mock[r["sample"]][canon(r["species"])] += int(r["reads"])

    tot = collections.Counter()
    for c in reads.values():
        for a, n in c.items():
            tot[kind[a]] += n
    all_reads = sum(tot.values())
    lines = [f"{name}: {len(reads)} samples, {len(asvs)} ASVs, {all_reads:,} reads; "
             + ", ".join(f"{k} {100 * v / all_reads:.1f}%" for k, v in tot.most_common()), ""]
    def clean_errors(calls):
        # per-sample species abundance error: clean ASVs; all reads vs the mock
        errs = []
        for smp in reads.values():
            t, p = collections.Counter(), collections.Counter()
            m = sum(v for a, v in smp.items() if kind[a] == "clean")
            for a, v in smp.items():
                if kind[a] != "clean":
                    continue
                t[truth[a]] += v / m
                rank, nm = calls[a]
                p[canon(nm) if rank == "species" else "(not species)"] += v / m
            errs.append(0.5 * sum(abs(t[k] - p[k]) for k in set(t) | set(p)))
        errm = []
        for smp_name, smp in reads.items():
            mt = mock[smp_name]
            mn, an = sum(mt.values()), sum(smp.values())
            p = collections.Counter()
            for a, v in smp.items():
                rank, nm = calls[a]
                p[canon(nm) if rank == "species" else "(not species)"] += v / an
            errm.append(0.5 * sum(abs(mt[k] / mn - p[k]) for k in set(mt) | set(p)))
        return f"  {100 * sum(errs) / len(errs):9.1f}%  {100 * sum(errm) / len(errm):6.1f}%"

    for which, title in (("clean", "Clean ASVs (an exact source copy; truth = the mock label)"),
                         ("noisy", f"Noisy ASVs (1-{NEAR} edits from sources of one species)")):
        if not any(kind[a] == which for a in asvs):
            continue
        lines.append(f"{title}, read-weighted:")
        lines.append(f"   {'arm':28s}" + "".join(f"{c:>17s}" for c in metrics)
                     + ("  abund. err  vs mock" if which == "clean" else ""))
        for label, calls in arms.items():
            c = collections.Counter()
            for smp in reads.values():
                for a, n in smp.items():
                    if kind[a] == which:
                        c[outcome(calls[a], truth[a])] += n
                        if calls[a][0] in ("species", "genus"):
                            c["genus calls"] += n
                            c["genus right"] += n * (call_genus(calls[a]) == truth[a].split()[0])
            n = sum(c[x] for x in cats)
            sp_made = c["right species"] + c["WRONG species"]
            vals = [c["right species"] / n, c["right species"] / sp_made if sp_made else float("nan"),
                    c["genus right"] / n, c["genus right"] / c["genus calls"] if c["genus calls"] else float("nan")]
            row = f"   {label:28s}" + "".join(f"{100 * v:16.1f}%" for v in vals)
            if which == "clean":
                row += clean_errors(calls)
            lines.append(row)
        lines.append("")
    lines.append("")
    for label, calls in arms.items():
        wrong = collections.Counter()
        for smp in reads.values():
            for a, n in smp.items():
                if kind[a] in ("clean", "noisy") and outcome(calls[a], truth[a]) in ("WRONG species", "WRONG genus"):
                    wrong[(truth[a], calls[a][1])] += n
        if wrong:
            lines.append(f"   {label} wrong: " + "; ".join(
                f"{t} -> {c} {100 * n / all_reads:.2f}%" for (t, c), n in wrong.most_common(6)))
    cur = {a for a in asvs if kind[a] == "curated"}
    if cur:
        lines.append("")
        lines.append("Curated sources (mock label is a record this database drops as mislabelled):")
        for label, calls in arms.items():
            c = collections.Counter()
            for smp in reads.values():
                for a, n in smp.items():
                    if a in cur:
                        c[f"{truth[a]} -> {canon(calls[a][1]) or 'unassigned'} ({calls[a][0] or '-'})"] += n
            lines.append(f"   {label:28s} " + "; ".join(f"{k} {100 * v / all_reads:.2f}%"
                                                        for k, v in c.most_common(4)))
    text = "\n".join(lines)
    open(P + ".score.txt", "w").write(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
