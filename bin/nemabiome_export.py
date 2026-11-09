#!/usr/bin/env python3
"""Turn the nemabiome web app's Excel export into backbone files.

    bin/nemabiome_export.py <database.xlsx> <out prefix>

The live database (https://ucvm-bio.shinyapps.io/Nematode_ITS2, embedded
at nemabiome.ca/its2-database) is only downloadable through the app, as a
zip holding database.xlsx and/or dada2.fasta. Export it with every column
shown, "sequence" included - the Excel file holds only the columns on
screen. Writes the same pair Zenodo v1.0.0 shipped, so the rest of the
build reads either:

  <prefix>_seqs_nr.fasta      >accn
  <prefix>_taxonomy_nr.txt    uid taxid accn superkingdom kingdom phylum
                              class order family genus species

(the app calls the first rank "domain"; renamed to match v1.0.0).
Reads the xlsx with the standard library only.
"""

from __future__ import annotations

import re
import sys
import zipfile
from xml.etree import ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def col_index(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1


def rows(path):
    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
        sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
    for row in sheet.iter(f"{{{NS['m']}}}row"):
        out = {}
        for c in row.findall("m:c", NS):
            t = c.get("t")
            if t == "inlineStr":
                v = "".join(x.text or "" for x in c.iter(f"{{{NS['m']}}}t"))
            else:
                e = c.find("m:v", NS)
                v = "" if e is None else e.text or ""
                if t == "s":
                    v = shared[int(v)]
            out[col_index(c.get("r"))] = v
        yield [out.get(i, "") for i in range(max(out) + 1)] if out else []


def main() -> int:
    xlsx, prefix = sys.argv[1:3]
    it = rows(xlsx)
    head = ["superkingdom" if h == "domain" else h for h in next(it)]
    need = ["accn", "uid", "taxid", "sequence"] + RANKS
    missing = [h for h in need if h not in head]
    if missing:
        sys.exit(f"export lacks columns {missing}: show every column in the app before exporting")
    ix = {h: head.index(h) for h in need}
    n = 0
    with open(prefix + "_seqs_nr.fasta", "w") as fa, open(prefix + "_taxonomy_nr.txt", "w") as tx:
        tx.write("\t".join(["uid", "taxid", "accn"] + RANKS) + "\n")
        for r in it:
            if not r:
                continue
            r = r + [""] * (len(head) - len(r))
            fa.write(f">{r[ix['accn']]}\n{r[ix['sequence']].upper()}\n")
            tx.write("\t".join(r[ix[h]] for h in ["uid", "taxid", "accn"] + RANKS) + "\n")
            n += 1
    print(f"{n} records", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
