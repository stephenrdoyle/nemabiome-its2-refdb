#!/usr/bin/env python3
"""GenBank ITS2 records the backbone lacks, in the backbone's own formats.

    bin/genbank_update.py <backbone taxonomy_nr.txt> <build dir> <since YYYY/MM/DD> <out dir> [accessions.txt]
    bin/genbank_update.py merge <backbone prefix> <backbone records.fasta.gz> <update dir> <out dir>

The backbone (nemabiome 20260317) was cut from NCBI nt by markerDB's HMM
and collapsed to unique sequences. This fetches what has come in since,
and what was never annotated as ITS2, so bin/its2.py and
bin/its2_partials.py can cut and grade it exactly as they did the
backbone. Three Entrez queries, nematodes only:

  new ITS2     ITS2 in the record, published since <since>
  new rDNA     5.8S, 28S or large subunit, published since <since>
  unlabelled   5.8S and 28S / large subunit in one record, never called
               ITS2 (any date): the ITS2 lies between them
  informal     ITS2 in the record (any date), named "Genus sp. X": the
               backbone holds no such record, and some lineages are
               known only so (the gorilla Necator of Hasegawa et al.)
  named        the accessions in [accessions.txt] (one per line, # notes),
               records found missing from the backbone by hand. Any
               informal name is kept for these

Kept: records whose accession is not in the backbone, 80-20,000 bases,
whose taxon sits at or below an NCBI species with a formal name, or is a
"Genus sp. X" under a nematode genus. NCBI ranks "Genus sp. X", cf.,
aff. and hybrids as species; the build gives a "sp." record its genus
only (species_status unnamed), and curation/lineages.tsv can carry the
name. cf., aff., hybrids, environmental and uncultured samples, and
"sp." names above genus ("Strongyloidea sp.") are left out. Lineage: NCBI Taxonomy's superkingdom (or domain) ...
species. The ITS2 region of each record, which bin/its2.py needs as the
backbone's HMM cut, is where it aligns to the build's records (minimap2,
the best hit by matches less twice the differences, extended by the
reference's unaligned ends): >= 80
aligned bases at >= 70% identity, else the record holds no ITS2 and is
left out.

merge writes the backbone and the update as one input for bin/its2.py
and bin/build.py: seqs.fasta, records.fasta.gz and taxonomy.txt (the
backbone's rows with a blank source, then the update's).

Writes to <out dir>: records.fasta.gz (full GenBank records, headers
"accession.version title"), seqs.fasta (the ITS2 region, for its2.py),
taxonomy.txt (backbone columns plus source), and log.txt (counts at each
step). Uses NCBI E-utilities at their documented rate; set NCBI_API_KEY
to go faster.
"""
from __future__ import annotations

import collections
import csv
import gzip
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build import read_fasta  # noqa: E402
from seqmatch import revcomp  # noqa: E402

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
KEY = os.environ.get("NCBI_API_KEY", "")
PAUSE = 0.12 if KEY else 0.4
SOURCE = "genbank_update"
RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
ITS2 = '("internal transcribed spacer 2"[All Fields] OR ITS2[All Fields] OR ITS-2[All Fields])'
MIN_LEN, MAX_LEN = 80, 20000
MIN_ALN, MIN_ID = 80, 0.70
INFORMAL = re.compile(r" sp\.| spp\.| cf\.| aff\.| x |[0-9]|^[A-Z][a-z]+ [A-Z]|uncultured|environmental")
NOT_SP = re.compile(r" cf\.| aff\.| x |uncultured|environmental")


def genus_sp(ranks):
    """An informal "Genus sp. X" name under a genus: kept at genus."""
    sp, genus = ranks["species"], ranks["genus"]
    return bool(genus) and re.match(re.escape(genus) + r" spp?\. ", sp + " ") is not None \
        and not NOT_SP.search(sp)


def call(util, **params):
    if KEY:
        params["api_key"] = KEY
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(5):
        try:
            with urllib.request.urlopen(urllib.request.Request(E + util, data=data), timeout=120) as r:
                out = r.read()
            time.sleep(PAUSE)
            return out
        except Exception:
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"{util} failed: {params.get('term', '')[:80]}")


def search(term):
    j = json.loads(call("esearch.fcgi", db="nuccore", term=term, retmax=100000, retmode="json"))
    return j["esearchresult"]["idlist"]


def summaries(uids):
    out = {}
    for k in range(0, len(uids), 400):
        d = json.loads(call("esummary.fcgi", db="nuccore", id=",".join(uids[k:k + 400]),
                            retmode="json"))["result"]
        for u in d["uids"]:
            r = d[u]
            out[u] = (r["accessionversion"], str(r["taxid"]), int(r["slen"]), r["title"])
    return out


def lineages(taxids, informal_ok=frozenset()):
    """taxid -> {rank: name} for the eight ranks, species from the species-rank
    ancestor (or the taxon itself); None where no species rank exists, or
    where it is an informal name other than "Genus sp. X" and the taxid is
    not in informal_ok."""
    out = {}
    for k in range(0, len(taxids), 200):
        root = ET.fromstring(call("efetch.fcgi", db="taxonomy", id=",".join(taxids[k:k + 200])))
        for t in root.findall("Taxon"):
            tid = t.findtext("TaxId")
            ranks = {}
            for a in t.findall("LineageEx/Taxon"):
                ranks[a.findtext("Rank")] = a.findtext("ScientificName")
            ranks[t.findtext("Rank")] = t.findtext("ScientificName")
            if "superkingdom" not in ranks and "domain" in ranks:
                ranks["superkingdom"] = ranks["domain"]
            lin = {r: ranks.get(r, "") for r in RANKS}
            out[tid] = (lin if lin["species"] and lin["phylum"] == "Nematoda"
                        and (tid in informal_ok or not INFORMAL.search(lin["species"])
                             or genus_sp(lin))
                        else None)
    return out


def fetch_fasta(accs):
    seqs = {}
    for k in range(0, len(accs), 200):
        text = call("efetch.fcgi", db="nuccore", id=",".join(accs[k:k + 200]), rettype="fasta",
                    retmode="text").decode()
        head = None
        for line in text.splitlines():
            if line.startswith(">"):
                head = line[1:].strip()
                seqs[head] = []
            elif line.strip() and head:
                seqs[head].append(line.strip().upper())
    return {h: "".join(v) for h, v in seqs.items()}


def its2_regions(records, refdb):
    """accession -> the record's ITS2 region, by its best minimap2 hit among the
    build's records, extended by the reference's unaligned ends."""
    with tempfile.TemporaryDirectory() as d:
        q = os.path.join(d, "q.fa")
        with open(q, "w") as fh:
            fh.writelines(f">{a}\n{s}\n" for a, s in records.items())
        paf = subprocess.run(["minimap2", "-c", "-t", "8", "-k", "11", "-w", "5", "-f", "1000000",
                              "--secondary=no", refdb, q], capture_output=True, text=True,
                             check=True).stdout
    best = {}
    for line in paf.splitlines():
        f = line.split("\t")
        qn, ql, qs, qe, st, tl, ts, te, m, al = (f[0], int(f[1]), int(f[2]), int(f[3]), f[4],
                                                 int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10]))
        if al < MIN_ALN or m / al < MIN_ID:
            continue
        # Matches less twice the differences: a long, loose hit (a backbone
        # record that runs into ITS1) must not beat an exact ITS2 one.
        sc = m - 2 * (al - m)
        if qn not in best or sc > best[qn][0]:
            best[qn] = (sc, ql, qs, qe, st, tl, ts, te)
    out = {}
    for a, (m, ql, qs, qe, st, tl, ts, te) in best.items():
        s = records[a]
        if st == "+":
            lo, hi = max(0, qs - ts), min(ql, qe + (tl - te))
            out[a] = s[lo:hi]
        else:
            lo, hi = max(0, qs - (tl - te)), min(ql, qe + ts)
            out[a] = revcomp(s[lo:hi])
    return out


def merge(prefix, records, update, outdir):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "seqs.fasta"), "w") as fh:
        for path in (prefix + "_seqs_nr.fasta", os.path.join(update, "seqs.fasta")):
            fh.write(open(path).read())
    with open(os.path.join(outdir, "records.fasta.gz"), "wb") as fh:
        for path in (records, os.path.join(update, "records.fasta.gz")):
            fh.write(open(path, "rb").read())
    bb = list(csv.DictReader(open(prefix + "_taxonomy_nr.txt"), delimiter="\t"))
    up = list(csv.DictReader(open(os.path.join(update, "taxonomy.txt")), delimiter="\t"))
    with open(os.path.join(outdir, "taxonomy.txt"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(bb[0]) + ["source"], delimiter="\t",
                           lineterminator="\n", extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(bb + up)
    print(f"{len(bb)} backbone + {len(up)} update records -> {outdir}")
    return 0


def main() -> int:
    if sys.argv[1] == "merge":
        return merge(*sys.argv[2:6])
    taxonomy, build, since, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    log = []
    backbone = {r["accn"].split(".")[0] for r in csv.DictReader(open(taxonomy), delimiter="\t")}
    nem, date = "Nematoda[Organism]", f"{since}:3000[PDAT]"
    queries = {
        "new ITS2": f"{nem} AND {ITS2} AND {date}",
        "new rDNA": f'{nem} AND (5.8S[All Fields] OR 28S[All Fields] OR "large subunit"[All Fields])'
                    f" NOT {ITS2} AND {date}",
        "unlabelled": f'{nem} AND 5.8S[All Fields] AND (28S[All Fields] OR "large subunit"[All Fields])'
                      f" NOT {ITS2}",
    }
    named = []
    if len(sys.argv) > 5:
        named = [ln.split()[0] for ln in open(sys.argv[5]) if ln.strip() and not ln.startswith("#")]
        queries["named"] = " OR ".join(f"{a}[Accession]" for a in named)
    # Last: a record any query above finds is credited to that one.
    queries["informal"] = f"{nem} AND {ITS2}"
    uids, why = [], {}
    for name, term in queries.items():
        ids = search(term)
        log.append(f"{name}: {len(ids)} records  ({term[:300]})")
        for u in ids:
            if u not in why:
                why[u] = name
                uids.append(u)
    summ = summaries(uids)
    keep = {u: s for u, s in summ.items() if s[0].split(".")[0] not in backbone}
    log.append(f"not in the backbone: {len(keep)} of {len(summ)}")
    keep = {u: s for u, s in keep.items() if MIN_LEN <= s[2] <= MAX_LEN}
    log.append(f"{MIN_LEN}-{MAX_LEN} bases: {len(keep)}")
    named_base = {a.split(".")[0] for a in named}
    lin = lineages(sorted({s[1] for s in keep.values()}),
                   frozenset(s[1] for s in keep.values() if s[0].split(".")[0] in named_base))
    def wanted(u, s):
        if s[0].split(".")[0] in named_base:
            return True
        formal = not INFORMAL.search(lin[s[1]]["species"])
        return (not formal if why[u] == "informal" else formal) or genus_sp(lin[s[1]])
    keep = {u: s for u, s in keep.items() if lin.get(s[1]) and wanted(u, s)}
    log.append(f"at or below a formally named nematode species, or Genus sp. X: {len(keep)}")
    fasta = fetch_fasta([s[0] for s in keep.values()])
    records = {h.split()[0].split(".")[0]: v for h, v in fasta.items()}
    titles = {h.split()[0].split(".")[0]: h for h in fasta}
    regions = its2_regions(records, os.path.join(build, "refdb.fasta"))
    log.append(f"holding an ITS2 region (>= {MIN_ALN} bases at >= {MIN_ID:.0%} to a build record): "
               f"{len(regions)}")
    by_why = collections.Counter(why[u] for u, s in keep.items() if s[0].split(".")[0] in regions)
    log.append("  by query: " + ", ".join(f"{k} {v}" for k, v in by_why.items()))
    with gzip.open(os.path.join(outdir, "records.fasta.gz"), "wt") as fh:
        for a in regions:
            fh.write(f">{titles[a]}\n{records[a]}\n")
    with open(os.path.join(outdir, "seqs.fasta"), "w") as fh:
        for u, s in keep.items():
            a = s[0].split(".")[0]
            if a in regions:
                fh.write(f">{s[0]}\n{regions[a]}\n")
    with open(os.path.join(outdir, "taxonomy.txt"), "w") as fh:
        fh.write("\t".join(["uid", "taxid", "accn"] + RANKS + ["source", "query"]) + "\n")
        for u, s in keep.items():
            if s[0].split(".")[0] in regions:
                fh.write("\t".join([u, s[1], s[0]] + [lin[s[1]][r] for r in RANKS]
                                   + [f"{SOURCE}_{since.replace('/', '')}", why[u]]) + "\n")
    text = "\n".join(log) + "\n"
    open(os.path.join(outdir, "log.txt"), "w").write(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
