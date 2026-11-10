#!/usr/bin/env python3
"""Current NCBI Taxonomy names for the backbone's taxids.

    bin/ncbi_names.py <backbone taxonomy_nr.txt> <out names.tsv>

markerDB named each record when the backbone was built; NCBI has renamed
taxa since (a genus moved, a species transferred, a synonym merged). This
fetches every taxid's current record and writes taxid, then NCBI's
current order, family, genus and species (empty where NCBI has no such
rank above the taxon), and the date fetched. bin/build.py --names
applies it: a record whose names differ takes NCBI's at those four ranks
(curation/overrides.tsv still has the last word). A merged taxid
resolves to the taxon it was merged into.
"""

from __future__ import annotations

import csv
import datetime
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
RANKS = ["order", "family", "genus", "species"]
BATCH = 200


def main() -> int:
    src, out = sys.argv[1:3]
    tids = sorted({r["taxid"] for r in csv.DictReader(open(src), delimiter="\t")}, key=int)
    cur = {}
    for i in range(0, len(tids), BATCH):
        data = urllib.parse.urlencode({"db": "taxonomy", "id": ",".join(tids[i:i + BATCH])}).encode()
        x = ET.fromstring(urllib.request.urlopen(EUTILS, data=data, timeout=120).read())
        time.sleep(0.4)
        for t in x.findall("Taxon"):
            lin = {t.findtext("Rank"): t.findtext("ScientificName")}
            for a in t.findall("LineageEx/Taxon"):
                lin.setdefault(a.findtext("Rank"), a.findtext("ScientificName"))
            names = [lin.get(r, "") for r in RANKS]
            for k in [t.findtext("TaxId")] + [a.text for a in t.findall("AkaTaxIds/TaxId")]:
                cur[k] = names
    missing = [t for t in tids if t not in cur]
    if missing:
        sys.exit(f"{len(missing)} taxids not returned by NCBI: {missing[:10]}")
    today = datetime.date.today().isoformat()
    with open(out, "w") as fh:
        fh.write("taxid\t" + "\t".join(RANKS) + "\tfetched\n")
        for t in tids:
            fh.write(t + "\t" + "\t".join(cur[t]) + f"\t{today}\n")
    print(f"{len(tids)} taxids -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
