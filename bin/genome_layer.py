#!/usr/bin/env python3
"""ITS2 from genome assemblies and long reads, as a layer for bin/build.py.

    bin/genome_layer.py <sources.tsv> <backbone taxonomy_nr.txt> <out layer.tsv>

<sources.tsv>: one scan per line - scan file (from bin/its2_scan.py),
kind (assembly or reads), species, source name (the assembly or run
accession). Lines starting # are ignored.

  assembly  every distinct ITS2 copy is kept - except that in an array of
            >= BIG_ARRAY copies a variant seen once may be an assembly
            error, so there it must be seen twice. A small array cannot
            tell an error from a real copy, so all are kept (the source
            name says which assembly).
  reads     a variant is kept if it is >= MIN_SHARE of the run's ITS2
            reads and seen in >= MIN_READS of them: HiFi errors are random,
            real rDNA copies recur (as for V4 in vespa-18s-refdb).

Each species' eight ranks come from the backbone (a record of the species,
else of its genus), else from NCBI Taxonomy. Writes accession (source:
contig:position, or source:variant n), source, label, copies, sequence -
cut at the same NC1 and 28S anchors as the backbone, so the records merge
with it directly.
"""

from __future__ import annotations

import collections
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
MIN_SHARE, MIN_READS, MIN_LEN = 0.02, 5, 80
BIG_ARRAY = 20
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def ncbi_lineage(species):
    q = urllib.parse.urlencode({"db": "taxonomy", "term": f'"{species}"[Scientific Name]',
                                "retmode": "json"})
    ids = json.loads(urllib.request.urlopen(EUTILS + "esearch.fcgi?" + q, timeout=60).read())
    time.sleep(0.4)
    tid = ids["esearchresult"]["idlist"][0]
    x = ET.fromstring(urllib.request.urlopen(EUTILS + "efetch.fcgi?db=taxonomy&id=" + tid,
                                             timeout=60).read())
    time.sleep(0.4)
    got = {"species": species}
    for t in x.iter("Taxon"):
        r, n = t.findtext("Rank"), t.findtext("ScientificName")
        if r in RANKS and r not in got:
            got[r] = n
        if r == "clade" and n == "Metazoa":
            got.setdefault("kingdom", n)
    got.setdefault("superkingdom", "Eukaryota")
    got.setdefault("kingdom", "Metazoa")
    return [got.get(r, "") for r in RANKS]


def main() -> int:
    src_path, tax_path, out = sys.argv[1:4]
    by_sp, by_genus = {}, {}
    for r in csv.DictReader(open(tax_path), delimiter="\t"):
        lab = [r[k] for k in RANKS]
        by_sp.setdefault(r["species"], lab)
        by_genus.setdefault(r["genus"], lab[:7])
    lineage = {}
    rows, n_in = [], collections.Counter()
    for line in open(src_path):
        if not line.strip() or line.startswith("#"):
            continue
        scan, kind, species, source = line.rstrip("\n").split("\t")[:4]
        if species not in lineage:
            if species in by_sp:
                lineage[species] = by_sp[species]
            elif species.split()[0] in by_genus:
                lineage[species] = by_genus[species.split()[0]] + [species]
            else:
                lineage[species] = ncbi_lineage(species)
        hits = [r for r in csv.DictReader(open(scan), delimiter="\t") if int(r["length"]) >= MIN_LEN]
        n_in[source] = len(hits)
        v = collections.Counter(r["its2"] for r in hits)
        if kind == "assembly":
            where = {}
            for r in hits:
                where.setdefault(r["its2"], f"{r['seq_id']}:{r['position']}")
            keep = [(f"{source}:{where[s]}", s, c) for s, c in v.most_common()
                    if c >= 2 or len(hits) < BIG_ARRAY]
        else:
            n = len(hits)
            keep = [(f"{source}:v{k + 1}", s, c) for k, (s, c) in enumerate(v.most_common())
                    if c >= MIN_READS and c >= MIN_SHARE * n]
        for acc, s, c in keep:
            rows.append((acc, "long_read_assemblies" if kind == "assembly" else "long_reads",
                         ";".join(lineage[species]), c, s))
    with open(out, "w") as fh:
        fh.write("accession\tsource\tlabel\tcopies\tsequence\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    per = collections.Counter(r[2].split(";")[-1] for r in rows)
    print(f"{len(rows)} layer records, {len(per)} species, from {len(n_in)} scans", file=sys.stderr)
    for sp, n in sorted(per.items()):
        print(f"  {n:3d}  {sp}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
