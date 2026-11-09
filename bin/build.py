#!/usr/bin/env python3
"""Collapse the graded ITS2 records into a release directory.

    bin/build.py <its2 prefix> <taxonomy_nr.txt> <source name> <outdir> [--version V] [--workers N]

    bin/build.py work/its2/its2 \\
        sources/nemabiome_its2_20260317/Nematode_ITS2_20260317_taxonomy_nr.txt \\
        nemabiome_20260317 work/build

<its2 prefix> names the outputs of bin/its2.py (<prefix>_cut.tsv) and
bin/its2_partials.py (<prefix>_refs.fasta, <prefix>_refs.tsv).

Writes the same set of files as vespa-18s-refdb, so amphid reads either:

  refdb.fasta                   one record per unique ITS2, NC1 end to 28S
  taxonomy.tsv                  merged taxonomy, deepest rank, members,
                                grade and coverage (below)
  members.tsv                   every source record and what became of it
  conflicts.tsv                 records whose members disagree at order or
                                above, and how each was decided
  primer_fit.tsv                NC1/NC2 edits of each record's best member
  resolution_nanopore.tsv       deepest rank read whole, with margin
  resolution_illumina_<R>x<R>.tsv  deepest rank through the two read
                                windows, with margin
  species_margins.tsv           edits to the nearest record of another
                                species, whole and per geometry
  anchor_lookalikes.tsv         adapter-anchor near-copies (pre-flight mask)
  manifest.json, summary.txt

MERGING. Identical ITS2s become one record. At each rank the record
takes the single informative name its members agree on; two names stop
it there. Species names that do not name one species do not vote:
hybrids ("A x B"), cf./aff., sensu lato, environmental samples.

CONFLICTS. Members disagreeing at order or above (class, order: the
backbone is all Nematoda) cannot all be right - identical ITS2 in
Rhabditida and Trichinellida is a mislabel. A label carried by at least
two members and two thirds of them wins (majority); otherwise every
member is merged and the record stops where they agree (unresolved),
listed for review. Below order, disagreement is left as biology: the
record stops at the rank its members share.

GRADE. A record's grade is its best member's (bin/its2_partials.py):
complete, near_complete, partial_90, partial_75, partial, unmeasured;
coverage likewise (1 for complete). amphid scores a read over the span
its reference covers, so a partial reference is scored on fewer bases -
the grade says how many fewer.

WINDOWS. NC1 and NC2 are 20 bases. A 2xR read pair, primers trimmed,
sees the first R-20 bases of the ITS2 and the last R-20 bases of the
amplicon - of which ~47 are the 28S stretch the references leave out, so
R-67 of the ITS2. Records that short or shorter are seen whole.
"""

from __future__ import annotations

import argparse
import collections
import csv
import datetime
import hashlib
import json
import os
import re
import subprocess
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "bin"))
from its2 import NC1, NC2  # noqa: E402
from seqmatch import ANCHORS, LOOKALIKE_MAX_EDITS, find_anchor_lookalikes  # noqa: E402

RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
CONFLICT_RANKS = RANKS[:RANKS.index("order") + 1]
GRADES = ["complete", "near_complete", "partial_90", "partial_75", "partial", "unmeasured"]
PRIMER_LEN = 20
STRETCH_28S = 47
READ_LENGTHS = [150, 250]
TOP = 30

NOT_A_SPECIES = re.compile(r" x | (cf|aff)\.? |sensu lato| s\.l\.|environmental sample")


def windows(read_len):
    return read_len - PRIMER_LEN, read_len - PRIMER_LEN - STRETCH_28S


def informative(rank, name):
    if not name or name == "NA":
        return False
    return not (rank == "species" and NOT_A_SPECIES.search(name + " "))


def merge(labels):
    """Ranks all labels agree on, blank from the first disagreement down."""
    out, stopped = [], False
    for i, rank in enumerate(RANKS):
        vals = {lab[i] for lab in labels if informative(rank, lab[i])}
        if stopped or len(vals) != 1:
            stopped = True
            out.append("")
        else:
            out.append(vals.pop())
    return out


def merge_records(records):
    """The common ancestor of records the read windows cannot tell apart.
    Unlike merge, a blank is not a placeholder here but a record that
    already stopped there - its members disagreed - so it stops the group."""
    out, stopped = [], False
    for i in range(len(RANKS)):
        vals = {r[i] for r in records}
        if stopped or len(vals) != 1 or "" in vals:
            stopped = True
            out.append("")
        else:
            out.append(vals.pop())
    return out


def deepest(ranks):
    for rank, name in reversed(list(zip(RANKS, ranks))):
        if name:
            return rank, name
    return "", ""


def read_fasta(path):
    seqs, cur = {}, None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(line.strip().upper())
    return {k: "".join(v) for k, v in seqs.items()}


def infix_distance(ref, read):
    """Edits to place all of `ref` somewhere in `read`, the read's bases
    either side free: how amphid scores a read over its reference's span.
    Myers' bit-parallel algorithm (vespa-18s-refdb bin/vespa.py's
    edit_distance) without the boundary bit, so the read may start
    anywhere, and the best column taken, so it may end anywhere."""
    m = len(ref)
    if not m:
        return 0
    peq = collections.defaultdict(int)
    for i, c in enumerate(ref):
        peq[c] |= 1 << i
    full, top = (1 << m) - 1, 1 << (m - 1)
    pv, mv, score = full, 0, m
    best = m
    for c in read:
        eq = peq[c]
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & full)
        mh = pv & xh
        if ph & top:
            score += 1
        elif mh & top:
            score -= 1
        ph = (ph << 1) & full
        mh = (mh << 1) & full
        pv = mh | (~(xv | ph) & full)
        mv = ph & xv
        best = min(best, score)
    return best


def span_distance(a, b):
    """Edits between two records over the span the shorter-reaching one
    covers: each placed inside the other, the better taken. A partial is
    compared only on the bases it has, as a read is scored only over its
    reference's span; for two complete records this is global distance
    less any edits at the very ends."""
    return min(infix_distance(a, b), infix_distance(b, a))


# A window placed in the other record's window, with this much slack for
# indels earlier in the ITS2 shifting where the window falls.
WINDOW_SLACK = 10


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def git_commit():
    try:
        rev = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", HERE, "status", "--porcelain", "--", "bin"],
                               capture_output=True, text=True).stdout.strip()
        return rev + ("-dirty" if dirty else "")
    except Exception:  # noqa: BLE001
        return "uncommitted"


# --- margins ---------------------------------------------------------------

_S = _G = None


def _init(seqs, geoms):
    global _S, _G
    _S, _G = seqs, geoms


def _margin(job):
    rid, cands = job
    s = _S[rid]
    whole = sorted((span_distance(s, _S[c]), c) for c in cands)
    out = [rid, str(whole[0][0]) if whole else ""]
    near = [whole[0][1] if whole else ""]
    k = WINDOW_SLACK
    for w1, w2 in _G:
        if len(s) <= w1 + w2:
            out.append(out[1])
            near.append(near[0])
            continue
        best = None
        for _, c in whole:
            t = _S[c]
            d = (min(infix_distance(s[:w1], t[:w1 + k]), infix_distance(t[:w1], s[:w1 + k]))
                 + min(infix_distance(s[-w2:], t[-w2 - k:]), infix_distance(t[-w2:], s[-w2 - k:])))
            if best is None or d < best[0]:
                best = (d, c)
        out.append("" if best is None else str(best[0]))
        near.append("" if best is None else best[1])
    return out + near


def margins(fasta, seqs, tax, geoms, workers):
    """Edits from each species-level record to the nearest record naming
    another species (or none): minimap2 proposes the TOP most similar,
    span_distance decides. -f 1000000 so shared minimizers are not
    dropped across thousands of near-identical strongylid records; a
    candidate must chain matches over half the SHORTER record, so a
    partial still finds the complete records it sits inside."""
    p = subprocess.Popen(["minimap2", "-t", str(workers), "-X", "-k", "11", "-w", "5",
                          "-f", "1000000", "-N", "2000", "-p", "0.5", fasta, fasta],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    score = collections.defaultdict(dict)
    for line in p.stdout:
        f = line.split("\t", 11)
        q, t, m = f[0], f[5], int(f[9])
        if m < 0.5 * min(int(f[1]), int(f[6])):
            continue
        for a, b in ((q, t), (t, q)):
            if tax[a]["deepest_rank"] == "species" and tax[b]["species"] != tax[a]["species"]:
                d = score[a]
                if m > d.get(b, -1):
                    d[b] = m
                    if len(d) > 8 * TOP:
                        for c in sorted(d, key=d.get)[:4 * TOP]:
                            del d[c]
    if p.wait():
        sys.exit("minimap2 failed")
    jobs = [(rid, sorted(score[rid], key=lambda c: -score[rid][c])[:TOP])
            for rid, t in tax.items() if t["deepest_rank"] == "species"]
    with Pool(workers, initializer=_init, initargs=(seqs, geoms)) as pool:
        return {r[0]: r for r in pool.map(_margin, jobs, chunksize=100)}


def lookalike_chunk(items):
    return find_anchor_lookalikes(items, ANCHORS, LOOKALIKE_MAX_EDITS)


# --- main ------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("its2")
    ap.add_argument("taxonomy")
    ap.add_argument("source")
    ap.add_argument("outdir")
    ap.add_argument("--version", default="dev")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    out = a.outdir
    os.makedirs(out, exist_ok=True)

    refs = read_fasta(a.its2 + "_refs.fasta")
    grade = {r["accession"]: r for r in csv.DictReader(open(a.its2 + "_refs.tsv"), delimiter="\t")}
    cut = {r["accession"]: r for r in csv.DictReader(open(a.its2 + "_cut.tsv"), delimiter="\t")}
    label = {}
    for r in csv.DictReader(open(a.taxonomy), delimiter="\t"):
        label[r["accn"].split(".")[0]] = (r["accn"], [r[k] for k in RANKS])

    # --- collapse
    groups = collections.defaultdict(list)
    for acc, s in refs.items():
        groups[s].append(acc)
    seqs, tax, members, conflicts = {}, {}, [], []
    for s, accs in groups.items():
        rid = "ITS2_" + hashlib.sha1(s.encode()).hexdigest()[:12]
        labs = {acc: label[acc][1] for acc in accs}
        role = {acc: "kept" for acc in accs}
        top = merge(labs.values())
        bad = next((rk for rk in CONFLICT_RANKS if not top[RANKS.index(rk)]
                    and len({lab[RANKS.index(rk)] for lab in labs.values()}) > 1), None)
        if bad:
            i = RANKS.index(bad)
            votes = collections.Counter(tuple(lab[:i + 1]) for lab in labs.values())
            (win, n), = votes.most_common(1)
            if n >= 2 and n >= 2 / 3 * len(accs):
                decision = "majority"
                for acc, lab in labs.items():
                    if tuple(lab[:i + 1]) != win:
                        role[acc] = "overruled"
                top = merge([labs[acc] for acc in accs if role[acc] == "kept"])
            else:
                decision = "unresolved"
            conflicts.append({"id": rid, "decision": decision, "rank": bad, "members": len(accs),
                              "labels": " | ".join(f"{';'.join(k)} x{v}" for k, v in votes.most_common()),
                              "kept_label": ";".join(x for x in top if x)})
        rank, _ = deepest(top)
        gs = sorted((GRADES.index(grade[acc]["grade"]) for acc in accs if role[acc] == "kept"))
        cov = max((1.0 if grade[acc]["grade"] == "complete" else float(grade[acc]["coverage"] or 0))
                  for acc in accs)
        seqs[rid] = s
        tax[rid] = dict(zip(RANKS, top), id=rid, deepest_rank=rank, members=len(accs),
                        grade=GRADES[gs[0]], coverage="" if GRADES[gs[0]] == "unmeasured" else f"{cov:.3f}",
                        sources=a.source)
        for acc in accs:
            g = grade[acc]
            members.append({"id": rid, "source": a.source, "accession": label[acc][0],
                            "label": ";".join(labs[acc]), "role": role[acc], "grade": g["grade"],
                            "coverage": g["coverage"], "trimmed_start": g["trimmed_start"],
                            "trimmed_end": g["trimmed_end"]})
    order = sorted(seqs)

    with open(os.path.join(out, "refdb.fasta"), "w") as fh:
        for rid in order:
            fh.write(f">{rid}\n{seqs[rid]}\n")
    cols = ["id"] + RANKS + ["deepest_rank", "members", "grade", "coverage", "sources"]
    with open(os.path.join(out, "taxonomy.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(tax[rid] for rid in order)
    for name, rows in (("members.tsv", sorted(members, key=lambda r: (r["id"], r["accession"]))),
                       ("conflicts.tsv", sorted(conflicts, key=lambda r: r["id"]))):
        with open(os.path.join(out, name), "w") as fh:
            fields = list(rows[0]) if rows else ["id"]
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)

    # --- primer fit: the best-fitting member that holds both sites
    fit = {}
    for m in members:
        acc = m["accession"].split(".")[0]
        e1, e2 = cut[acc]["nc1_edits"], cut[acc]["nc2_edits"]
        known = (e1 != "") + (e2 != "")
        key = (-known, (int(e1) if e1 else 0) + (int(e2) if e2 else 0))
        if m["id"] not in fit or key < fit[m["id"]][0]:
            fit[m["id"]] = (key, e1, e2, m["accession"])
    fit_class = {}
    with open(os.path.join(out, "primer_fit.tsv"), "w") as fh:
        fh.write("id\tfit\tnc1_edits\tnc2_edits\tbest_member\n")
        for rid in order:
            _, e1, e2, acc = fit[rid]
            if e1 == "" and e2 == "":
                c = "unknown"
            elif e1 == "" or e2 == "":
                c = "one_site"
            elif int(e1) == int(e2) == 0:
                c = "exact"
            elif max(int(e1), int(e2)) <= 2:
                c = "near"
            else:
                c = "poor"
            fit_class[rid] = c
            fh.write(f"{rid}\t{c}\t{e1}\t{e2}\t{acc}\n")

    # --- margins
    geoms = [windows(r) for r in READ_LENGTHS]
    marg = margins(os.path.join(out, "refdb.fasta"), seqs, tax, geoms, a.workers)
    with open(os.path.join(out, "species_margins.tsv"), "w") as fh:
        fh.write("\t".join(["id", "margin_whole"] + [f"margin_{w1}x{w2}" for w1, w2 in geoms]
                           + ["nearest_whole"] + [f"nearest_{w1}x{w2}" for w1, w2 in geoms]) + "\n")
        for rid in sorted(marg):
            fh.write("\t".join(marg[rid]) + "\n")

    # --- resolution
    # A species-level record at margin 0 sits inside another species'
    # record over its own span: a read from it scores the same against
    # both, so it resolves only to the rank the two share ("tied").
    ranks_order = {r: i for i, r in enumerate(RANKS)}
    ng = len(geoms)

    def lab(i):
        return [tax[i][k] for k in RANKS]

    def tied(rid, col):
        m = marg.get(rid)
        if m and m[1 + col] == "0" and m[2 + ng + col]:
            return merge_records([lab(rid), lab(m[2 + ng + col])]), True
        return lab(rid), False

    nano, nano_tied, whole_lab = collections.Counter(), 0, {}
    with open(os.path.join(out, "resolution_nanopore.tsv"), "w") as fh:
        fh.write("id\tdeepest_rank\tname\tgrade\tmargin\n")
        for rid in order:
            whole_lab[rid], t = tied(rid, 0)
            nano_tied += t
            r, n = deepest(whole_lab[rid])
            nano[r] += 1
            mg = marg[rid][1] if rid in marg and r == "species" else ""
            fh.write(f"{rid}\t{r}\t{n}\t{tax[rid]['grade']}\t{mg}\n")
    geometries = []
    for gi, (read_len, (w1, w2)) in enumerate(zip(READ_LENGTHS, geoms)):
        key = {rid: (s if len(s) <= w1 + w2 else s[:w1] + "|" + s[-w2:]) for rid, s in seqs.items()}
        by = collections.defaultdict(list)
        for rid, k in key.items():
            by[k].append(rid)
        gid = {k: "W_" + hashlib.sha1(k.encode()).hexdigest()[:12] for k in by}
        table = f"resolution_illumina_{read_len}x{read_len}.tsv"
        ill, ill_tied = collections.Counter(), 0
        with open(os.path.join(out, table), "w") as fh:
            fh.write("id\twindow_group\trecords_in_group\tdeepest_rank\tname\tgrade\tmargin\n")
            for rid in order:
                ids = by[key[rid]]
                win, t = tied(rid, 1 + gi)
                ill_tied += t
                # Never deeper than read whole: fold in the whole-read call.
                labs = [lab(i) for i in ids] + [win, whole_lab[rid]]
                r, n = deepest(merge_records(labs))
                whole = deepest(whole_lab[rid])[0]
                if ranks_order.get(r, -1) > ranks_order.get(whole, -1):
                    sys.exit(f"{rid} resolves deeper through the windows than whole")
                ill[r] += 1
                mg = marg[rid][2 + gi] if rid in marg and r == "species" else ""
                fh.write(f"{rid}\t{gid[key[rid]]}\t{len(ids)}\t{r}\t{n}\t{tax[rid]['grade']}\t{mg}\n")
        geometries.append({"read_length": read_len, "windows": [w1, w2], "file": table,
                           "window_groups": len(by), "deepest_rank": dict(ill),
                           "tied": ill_tied})

    # --- anchor look-alikes
    items = [(rid, seqs[rid]) for rid in order]
    chunks = [items[i:i + 500] for i in range(0, len(items), 500)]
    with Pool(a.workers) as pool:
        look = [r for part in pool.map(lookalike_chunk, chunks) for r in part]
    with open(os.path.join(out, "anchor_lookalikes.tsv"), "w") as fh:
        fh.write("seq_id\tanchor\tstrand\tposition\tedits\tcontext\n")
        for r in look:
            fh.write(f"{r['seq_id']}\t{r['anchor']}\t{r['strand']}\t{r['position']}\t{r['edits']}\t{r['context']}\n")

    # --- manifest and summary
    gcount = collections.Counter(tax[rid]["grade"] for rid in order)
    roles = collections.Counter(m["role"] for m in members)
    dec = collections.Counter(c["decision"] for c in conflicts)
    pf = collections.Counter(fit_class.values())
    named = {t["species"] for t in tax.values() if t["species"]}
    manifest = {
        "name": f"nemits2_{a.version}", "version": a.version,
        "built": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "code_commit": git_commit(), "marker": "ITS2 (nematodes)",
        "primers": {"forward": {"name": "NC1", "sequence": NC1},
                    "reverse": {"name": "NC2", "sequence": NC2},
                    "reference": "Gasser et al. 1993, Nucleic Acids Res 21:2525"},
        "insert": "from the 3' end of NC1 to the start of the 28S stretch before NC2; "
                  "primers and the ~47-base 28S stretch not included",
        "score_over": "the reference's own span (a read's 28S stretch counts neither for nor against)",
        "illumina_geometries": [{k: g[k] for k in ("read_length", "windows", "file", "window_groups")}
                                for g in geometries],
        "anchors": ANCHORS, "lookalike_max_edits": LOOKALIKE_MAX_EDITS,
        "source": a.source, "source_taxonomy_sha256": sha256(a.taxonomy),
        "counts": {"records": len(order), "members": len(members), "named_species": len(named),
                   "grade": dict(gcount), "deepest_rank_nanopore": dict(nano),
                   "member_roles": dict(roles), "conflicts": dict(dec), "primer_fit": dict(pf)},
        "files": {},
    }
    for f in sorted(os.listdir(out)):
        if f not in ("manifest.json", "summary.txt"):
            manifest["files"][f] = sha256(os.path.join(out, f))

    def table(c, keys, n):
        return [f"    {k or '(none)':14s} {c[k]:7,d}  {100 * c[k] / n:5.1f}%" for k in keys]
    keys = list(reversed(RANKS)) + [""]
    lines = [f"{manifest['name']}  built {manifest['built']}  code {manifest['code_commit']}",
             f"primers NC1 {NC1} / NC2 {NC2}", f"source {a.source}",
             f"records {len(order):,} from {len(members):,} source records; "
             f"named species {len(named):,}",
             "members: " + ", ".join(f"{k} {v:,}" for k, v in sorted(roles.items())),
             "conflicts (order or above): " + (", ".join(f"{k} {v}" for k, v in sorted(dec.items())) or "none"),
             "", "Grade (best member):"] + table(gcount, GRADES, len(order))
    lines += ["", f"Deepest rank - nanopore (whole ITS2; {nano_tied:,} species-level records "
                  f"tied with another species' over their span):"] + table(nano, keys, len(order))
    for g in geometries:
        lines += [f"  Illumina {g['read_length']}x{g['read_length']} (windows {g['windows'][0]} + "
                  f"{g['windows'][1]}, {g['window_groups']:,} window groups; {g['tied']:,} tied):"]
        lines += table(collections.Counter(g["deepest_rank"]), keys, len(order))
    lines += ["", "Primer fit (best member holding the sites):"] + table(
        pf, ["exact", "near", "poor", "one_site", "unknown"], len(order))
    mw = [int(r[1]) for r in marg.values() if r[1] != ""]
    lines += ["", f"Species margins over the shorter record's span, read whole "
                  f"({len(marg):,} species-level records):"]
    for lo, hi, lab in ((0, 0, "0 (inside another species' record: tied)"), (1, 1, "1 edit"),
                        (2, 2, "2 edits"), (3, 5, "3-5 edits"), (6, 10, "6-10 edits"),
                        (11, 10 ** 6, "11+ edits")):
        n = sum(lo <= x <= hi for x in mw)
        lines.append(f"    {lab:44s} {n:7,d}")
    lines.append(f"    {'no other species within reach':44s} {len(marg) - len(mw):7,d}")
    text = "\n".join(lines) + "\n"
    open(os.path.join(out, "summary.txt"), "w").write(text)
    manifest["files"]["summary.txt"] = sha256(os.path.join(out, "summary.txt"))
    with open(os.path.join(out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
        fh.write("\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
