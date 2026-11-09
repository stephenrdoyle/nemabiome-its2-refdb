#!/usr/bin/env python3
"""Fetch the full GenBank record behind every backbone ITS2 sequence.

    bin/fetch_backbone_records.py <taxonomy_nr.txt> <out.fasta.gz>

The nemabiome backbone holds ITS2 only, cut by an HMM whose boundaries
drift from record to record. Re-cutting at fixed anchors (bin/its2.py)
needs the 5.8S and 28S around it, which only the source records have.
Fetches each accn in batches of 200 and writes them gzipped, headers
"accession.version title". Resumable: accessions already in the output
are skipped and new ones appended.

Uses NCBI E-utilities at their documented rate (3 requests a second
without an API key); set NCBI_API_KEY to go faster.
"""

from __future__ import annotations

import csv
import gzip
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
KEY = os.environ.get("NCBI_API_KEY", "")
PAUSE = 0.12 if KEY else 0.4
BATCH = 200


def efetch(ids: list[str]) -> str:
    params = {"db": "nuccore", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"}
    if KEY:
        params["api_key"] = KEY
    data = urllib.parse.urlencode(params).encode()
    for attempt in range(8):
        try:
            with urllib.request.urlopen(EUTILS, data=data, timeout=180) as r:
                body = r.read().decode()
            time.sleep(PAUSE)
            return body
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(15 * (attempt + 1))
    raise RuntimeError(f"NCBI efetch failed eight times: {err!r}")


def main() -> int:
    tax, out = sys.argv[1:3]
    want = list(dict.fromkeys(r["accn"] for r in csv.DictReader(open(tax), delimiter="\t")))
    have = set()
    if os.path.exists(out):
        with gzip.open(out, "rt") as fh:
            have = {line[1:].split()[0].split(".")[0] for line in fh if line.startswith(">")}
    todo = [a for a in want if a.split(".")[0] not in have]
    print(f"{len(want)} accessions, {len(have)} already fetched, {len(todo)} to go", file=sys.stderr)
    with gzip.open(out, "at") as fh:
        for i in range(0, len(todo), BATCH):
            fh.write(efetch(todo[i:i + BATCH]).strip() + "\n")
            fh.flush()
            print(f"  {min(i + BATCH, len(todo))}/{len(todo)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
