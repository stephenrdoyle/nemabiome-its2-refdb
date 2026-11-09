#!/usr/bin/env python3
"""Find the ITS2 in a genome assembly or a set of raw reads.

    bin/its2_scan.py panel <its2_cut.tsv> <records.fasta.gz> <taxonomy_nr.txt> <panel.fa>
    bin/its2_scan.py scan <panel.fa> <source> <out.tsv> [--reads N] [--threads T]

<source>: a path or URL (http/https/ftp; EBI's FTP fallback as in
vespa-18s-refdb's read scans), FASTA or FASTQ, plain or gzipped.

Why not a primer scan as for V4: NC1 is 4-6 bases off in many orders
outside the Strongylida, and no 5.8S k-mer is shared by all orders
(GAAGAACGCAGC is absent from every Trichinellida record). So the rDNA
is located by alignment, then cut by the anchors of bin/its2.py:

  panel  up to PER_GENUS source records per genus whose ITS2 bin/its2.py
         anchored at both ends - 5.8S, ITS2 and 28S stretch, from every
         genus in the backbone
  scan   assembly contigs are cut into WINDOW-base windows (overlapping
         by OVERLAP, so every amplicon sits whole in one); reads are
         taken whole. Each is mapped to the panel with minimap2 and only
         the hits kept; the hit region (with FLANK either side), turned
         to the panel's strand, is cut: the start at the 3' end of NC1
         (<= 3 edits, else <= 6), the end at the 28S stretch (<= 12
         edits) downstream, NC2 noted if present.

Writes one row per ITS2 found: source, sequence id, offset, strand, NC1 /
28S / NC2 edits, length, ITS2. Copies found twice in overlapping
windows are written once.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import io
import os
import re
import subprocess
import sys
import threading
import time
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from its2 import NC1, NC2, S28, read_fasta  # noqa: E402
from seqmatch import find_anchor_cost, revcomp  # noqa: E402

PER_GENUS = 2
WINDOW, OVERLAP = 20000, 3000
FLANK = 400
MAX_ITS2 = 1500
MIN_ITS2 = 80
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


# --- panel -------------------------------------------------------------------

def build_panel(cut_tsv, records, taxonomy, out):
    tax = {r["accn"].split(".")[0]: r for r in csv.DictReader(open(taxonomy), delimiter="\t")}
    rec = read_fasta(records)
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(cut_tsv), delimiter="\t"):
        a = r["accession"]
        if r["start"] in ("primer", "loose") and r["end"] in ("primer", "loose") and a in tax:
            by[tax[a]["genus"]].append((r["nc2_edits"] == "", -int(r["length"]), a))
    n = 0
    with open(out, "w") as fh:
        for genus, accs in sorted(by.items()):
            # Prefer records that run on to NC2 (they carry the whole 28S
            # stretch), then longer ITS2s; keep only 5.8S tail to 28S head,
            # so 18S or ITS1 in a source record does not draw hits.
            for _, _, a in sorted(accs)[:PER_GENUS]:
                x = None
                for y in (rec[a], revcomp(rec[a])):
                    h1 = find_anchor_cost(y, NC1, 6)
                    h2 = h1 and find_anchor_cost(y[h1[1]:h1[1] + MAX_ITS2], S28, 12)
                    if h1 and h2:
                        x = y[max(0, h1[0] - 150):h1[1] + h2[1] + 80]
                        break
                if x:
                    fh.write(f">{a}_{genus}\n{x}\n")
                    n += 1
    print(f"{n} panel records from {len(by)} genera", file=sys.stderr)


# --- input -------------------------------------------------------------------

def open_source(src):
    if not src.startswith(("ftp.", "http://", "https://", "ftp://")):
        return open(src, "rb")
    path = src.split("://", 1)[-1]
    err = None
    for attempt, url in enumerate(["https://" + path, "https://" + path,
                                   "ftp://" + path, "ftp://" + path]):
        try:
            return urllib.request.urlopen(url, timeout=300)
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(10 * (attempt + 1))
    raise err


def records_of(src):
    """(name, sequence) from FASTA or FASTQ, plain or gzipped."""
    raw = open_source(src)
    head = raw.read(2)
    rest = io.BufferedReader(_Prepend(head, raw))
    fh = io.TextIOWrapper(gzip.GzipFile(fileobj=rest) if head == b"\x1f\x8b" else rest,
                          "ascii", errors="replace")
    first = fh.readline()
    if first.startswith("@"):
        line = first
        while line:
            name = line[1:].split()[0]
            seq = fh.readline().strip().upper()
            fh.readline()
            fh.readline()
            yield name, seq
            line = fh.readline()
    else:
        name, buf = first[1:].split()[0], []
        for line in fh:
            if line.startswith(">"):
                yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
        yield name, "".join(buf)


class _Prepend(io.RawIOBase):
    def __init__(self, head, raw):
        self.head, self.raw = head, raw

    def readable(self):
        return True

    def readinto(self, b):
        if self.head:
            n = min(len(b), len(self.head))
            b[:n], self.head = self.head[:n], self.head[n:]
            return n
        data = self.raw.read(len(b))
        b[:len(data)] = data
        return len(data)


# --- scan --------------------------------------------------------------------

def cut(x, loose_ok=True):
    """ITS2 in x (panel strand): (start, end, nc1, s28, nc2) or None.

    The matcher reads N as any base, so an assembly gap matches every
    anchor perfectly (Dracunculus: ten 50-base "ITS2"s of pure N). Gaps
    are turned into a character that matches nothing, and an ITS2 with
    any N in it is not reported.

    loose_ok: whether an NC1 site at 4-6 edits may be taken. A read that
    starts inside the ITS2 has no NC1 at all, and the best 6-edit
    look-alike inside the ITS2 then makes a truncated "ITS2" (11% of
    Bunostomum's reads, all 50 bases). The caller allows it only when
    the alignment starts at or before NC1 on the panel record, i.e. the
    5.8S end is really there. MIN_ITS2 catches what is left."""
    x = x.replace("N", "#")
    hit = find_anchor_cost(x, NC1, 3) or (loose_ok and find_anchor_cost(x, NC1, 6))
    if not hit:
        return None
    s = hit[1]
    end = find_anchor_cost(x[s + 50:s + MAX_ITS2], S28, 12)
    if not end:
        return None
    e = s + 50 + end[0]
    if "#" in x[s:e] or e - s < MIN_ITS2:
        return None
    n2 = find_anchor_cost(x[e:e + 90], revcomp(NC2), 6)
    return s, e, hit[2], end[2], (n2[2] if n2 else "")


def scan(panel, src, out, reads, threads):
    # Low chaining thresholds (-m/-s 20, -n 2): a lineage with no close
    # panel relative shares only the 5.8S tail and 28S head with it, and
    # minimap2's defaults dropped Spirocerca lupi's rDNA entirely. Spurious
    # hits cost nothing: the anchors decide.
    p = subprocess.Popen(["minimap2", "-a", "--sam-hit-only", "--secondary=no", "-k", "11",
                          "-w", "5", "-m", "20", "-s", "20", "-n", "2",
                          "-t", str(threads), panel, "-"],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, text=True)
    fed = collections.Counter()

    def feed():
        try:
            for i, (name, seq) in enumerate(records_of(src)):
                if reads:
                    if i >= reads:
                        break
                    p.stdin.write(f">{name}\n{seq}\n")
                    fed["reads"] += 1
                    continue
                fed["bases"] += len(seq)
                step = WINDOW - OVERLAP
                for off in range(0, max(1, len(seq) - OVERLAP), step):
                    p.stdin.write(f">{name}@{off}\n{seq[off:off + WINDOW]}\n")
        finally:
            p.stdin.close()
    t = threading.Thread(target=feed)
    t.start()
    nc1_end = {name: (find_anchor_cost(s, NC1, 6) or (0, 0, 0))[1]
               for name, s in read_fasta(panel).items()}
    seen, rows = set(), []
    for line in p.stdout:
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 4 or flag & 256:
            continue
        name, seq = f[0], f[9]
        if seq == "*":
            continue
        ops = CIGAR.findall(f[5])
        qs = int(ops[0][0]) if ops[0][1] in "SH" else 0
        alen = sum(int(n) for n, o in ops if o in "MIS=X") - qs - (
            int(ops[-1][0]) if ops[-1][1] in "SH" else 0)
        lo, hi = max(0, qs - FLANK), qs + alen + FLANK
        x = seq[lo:hi]
        c = cut(x, loose_ok=int(f[3]) - 1 <= nc1_end.get(f[2], 0))
        if not c:
            continue
        s, e, e1, e2, e3 = c
        contig, _, woff = name.partition("@")
        strand = "-" if flag & 16 else "+"
        # Position on the original contig/read, panel strand flipped back.
        L = len(seq)
        pos = (lo + s) if strand == "+" else (L - (lo + e))
        pos += int(woff or 0)
        key = (contig, pos)
        if key in seen:
            continue
        seen.add(key)
        rows.append((os.path.basename(src), contig, pos, strand, e1, e2, e3, e - s, x[s:e]))
    t.join()
    p.wait()
    with open(out, "w") as fh:
        fh.write("source\tseq_id\tposition\tstrand\tnc1_edits\ts28_edits\tnc2_edits\tlength\tits2\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    what = f"{fed['reads']:,} reads" if reads else f"{fed['bases'] / 1e6:.1f} Mb"
    print(f"{os.path.basename(src)}: {what}, {len(rows)} ITS2", file=sys.stderr)


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("panel")
    a.add_argument("cut_tsv")
    a.add_argument("records")
    a.add_argument("taxonomy")
    a.add_argument("out")
    b = sub.add_parser("scan")
    b.add_argument("panel")
    b.add_argument("source")
    b.add_argument("out")
    b.add_argument("--reads", type=int, default=0)
    b.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    if args.cmd == "panel":
        build_panel(args.cut_tsv, args.records, args.taxonomy, args.out)
    else:
        scan(args.panel, args.source, args.out, args.reads, args.threads)
    return 0


if __name__ == "__main__":
    sys.exit(main())
