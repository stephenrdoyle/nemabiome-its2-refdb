# nemabiome-its2-refdb

Builds the nematode ITS2 reference database for NC1/NC2 "nemabiome"
metabarcoding, read on Illumina (2x150, 2x250) and nanopore: every
nematode in the nemabiome ITS2 database, re-cut to the span the reads
cover, graded by how much of it each record holds, with taxonomy,
resolution tables and the adapter look-alike mask. It is the ITS2
counterpart of [vespa-18s-refdb](https://github.com/stephenrdoyle/vespa-18s-refdb),
writes the same release format, and is consumed by
[amphid](https://github.com/stephenrdoyle/nf-amphid-pipeline).

**Status: first scratch build. No release has been cut.** See
[docs/DESIGN.md](docs/DESIGN.md) for the design and what the first build
shows, and [docs/feedback/its2_notes_for_amphid.md](docs/feedback/its2_notes_for_amphid.md)
for what amphid needs to change to use it.

## What a record is

NC1 sits at the 3' end of 5.8S, NC2 in 28S. A record runs from the end
of NC1 to the start of the ~47-base 28S stretch before NC2 - the ITS2,
cut at fixed anchors rather than where the source database's HMM put
its boundaries. Reads carry the 28S stretch too; amphid scores a read
over its reference's span only, so it counts neither for nor against.

Most GenBank ITS2 records do not run primer to primer. Each record's
`grade` says how much of the ITS2 it holds, judged against its closest
record anchored at both ends: `complete`, `near_complete` (>= 98%),
`partial_90`, `partial_75`, `partial`, or `unmeasured` (no relative at
>= 80% identity).

## The first build (scratch, 2026-11-09)

12,988 records from 15,108 nemabiome sequences (version 20260317), 425
genera, 1,837 species with a record of their own.

| Grade | Records |
|---|---:|
| complete | 4,045 |
| near_complete | 3,709 |
| partial_90 | 1,302 |
| partial_75 | 971 |
| partial | 994 |
| unmeasured | 1,967 |

How far a record can be named, best case (an exact read, no sequencing
error). A record whose span sits exactly inside another species' record
is tied with it and stops at the rank they share.

| Read | Species | Genus | Family or above |
|---|---:|---:|---:|
| Nanopore (whole ITS2) | 93% | 6% | 1% |
| Illumina 150x150 | 92% | 7% | 1% |
| Illumina 250x250 | 93% | 6% | 1% |

Most species-level records are close to another species: 1,040 at one
edit, 859 at two, over a median ITS2 of 240 bases (`species_margins.tsv`).
At that length amphid's 99.7% species limit allows no mismatch at all.

## Build

```bash
S=sources/nemabiome_its2_20260317/Nematode_ITS2_20260317
python3 bin/nemabiome_export.py sources/nemabiome_its2_20260317/raw/database.xlsx $S
python3 bin/fetch_backbone_records.py ${S}_taxonomy_nr.txt sources/genbank_backbone/records.fasta.gz
python3 bin/its2.py ${S}_seqs_nr.fasta sources/genbank_backbone/records.fasta.gz work/its2/its2_cut
python3 bin/its2_partials.py work/its2/its2_cut work/its2/its2_refs 8
python3 bin/build.py work/its2/its2 ${S}_taxonomy_nr.txt nemabiome_20260317 work/build --workers 8
```

Needs Python 3.9+ and minimap2. `sources/README.md` says where each
download comes from.

## Release files

| file | holds |
|---|---|
| `refdb.fasta` | one record per unique ITS2 |
| `taxonomy.tsv` | superkingdom ... species, deepest rank, members, grade, coverage |
| `members.tsv` | every source record, its grade, and bases trimmed at each end |
| `conflicts.tsv` | records whose members disagree at order or above |
| `primer_fit.tsv` | NC1/NC2 edits of each record's best member, where its source record holds the sites |
| `resolution_nanopore.tsv`, `resolution_illumina_<R>x<R>.tsv` | deepest rank read whole or through the two read windows, grade, margin |
| `species_margins.tsv` | edits to the nearest record of another species, over the shorter record's span |
| `anchor_lookalikes.tsv` | adapter-anchor near-copies, for amphid's pre-flight mask |
| `manifest.json`, `summary.txt` | versions, primers, windows, counts, checksums |

## Citation

The backbone is the nemabiome ITS2 database (Workentine et al. 2020,
BMC Genetics 21:74; https://www.nemabiome.ca), CC BY 4.0. Primers:
Gasser et al. 1993, Nucleic Acids Res 21:2525; method: Avramenko et al.
2015, PLoS ONE 10:e0143559.

MIT licence for the code (LICENSE).
