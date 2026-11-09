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

12,973 records from 15,108 nemabiome sequences (version 20260317; 22
dropped as mislabelled, `curation/overrides.tsv`), 424 genera, 1,837
species with a record of their own.

<!-- figures:start -->
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/grade_dark.svg">
  <img alt="Bar chart: records per nematode group in nemits2_dev, by how much of the ITS2 they hold: I · animal parasites 401 (complete 24, near complete 176, partial 70, unmeasured 131); I · plant parasites 168 (complete 0, near complete 0, partial 0, unmeasured 168); I · free-living 16 (complete 2, near complete 2, partial 0, unmeasured 12); II · plant parasites 56 (complete 0, near complete 0, partial 0, unmeasured 56); II · free-living 19 (complete 1, near complete 0, partial 0, unmeasured 18); III · animal parasites 2,395 (complete 752, near complete 670, partial 759, unmeasured 214); IV · animal parasites 856 (complete 155, near complete 242, partial 153, unmeasured 306); IV · plant parasites 4,615 (complete 2,178, near complete 910, partial 946, unmeasured 581); IV · free-living 42 (complete 13, near complete 9, partial 6, unmeasured 14); V · animal parasites 4,063 (complete 771, near complete 1,633, partial 1,278, unmeasured 381); V · free-living 326 (complete 143, near complete 62, partial 51, unmeasured 70); Other Chromadoria · free-living 15 (complete 0, near complete 0, partial 0, unmeasured 15); Family not placed 1 (complete 1, near complete 0, partial 0, unmeasured 0)." src="docs/figures/grade_light.svg" width="760">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/taxa_dark.svg">
  <img alt="Bar chart: records per nematode group in nemits2_dev, and how many resolve to species read whole: I · animal parasites 401 (381 to species); I · plant parasites 168 (168 to species); I · free-living 16 (16 to species); II · plant parasites 56 (56 to species); II · free-living 19 (19 to species); III · animal parasites 2,395 (2,139 to species); IV · animal parasites 856 (806 to species); IV · plant parasites 4,615 (4,289 to species); IV · free-living 42 (42 to species); V · animal parasites 4,063 (3,901 to species); V · free-living 326 (314 to species); Other Chromadoria · free-living 15 (14 to species); Family not placed 1 (0 to species)." src="docs/figures/taxa_light.svg" width="760">
</picture>

How many species nemits2_dev can name. A species named in the source records but with no record of its own shares its ITS2 with another species (or sits inside another species' record over its span), so a read from it is named only to genus or above. Names are counted as written.

| Records | Species named in sources | Species with their own record (read whole) | Species named through Illumina 150x150 | Species named through Illumina 250x250 | Genera |
|---|---:|---:|---:|---:|---:|
| All records | 1,884 | 1,837 | 1,830 | 1,837 | 424 |
| **I · animal parasites** | 33 | 32 | 32 | 32 | 6 |
| I · plant parasites | 65 | 65 | 64 | 65 | 3 |
| I · free-living | 10 | 10 | 10 | 10 | 9 |
| II · plant parasites | 18 | 18 | 18 | 18 | 2 |
| II · free-living | 9 | 9 | 9 | 9 | 6 |
| **III · animal parasites** | 295 | 283 | 278 | 283 | 97 |
| **IV · animal parasites** | 152 | 147 | 147 | 147 | 13 |
| IV · plant parasites | 664 | 645 | 645 | 645 | 95 |
| IV · free-living | 21 | 21 | 21 | 21 | 11 |
| **V · animal parasites** | 510 | 506 | 505 | 506 | 156 |
| V · free-living | 104 | 100 | 100 | 100 | 24 |
| Other Chromadoria · free-living | 1 | 1 | 1 | 1 | 2 |
| Family not placed | 2 | 0 | 0 | 0 | 0 |

How far a record can be named in nemits2_dev, best case (an exact read, no sequencing error): the share of records resolving to species, and to genus at best.

| Records | Nanopore (whole ITS2): species | genus | Illumina 150x150: species | genus | Illumina 250x250: species | genus |
|---|---:|---:|---:|---:|---:|---:|
| All records (12,973) | 94% | 6% | 92% | 7% | 94% | 6% |
| **I · animal parasites** (401) | 95% | 5% | 90% | 10% | 94% | 6% |
| I · plant parasites (168) | 100% | 0% | 99% | 1% | 100% | 0% |
| I · free-living (16) | 100% | 0% | 100% | 0% | 100% | 0% |
| II · plant parasites (56) | 100% | 0% | 100% | 0% | 100% | 0% |
| II · free-living (19) | 100% | 0% | 100% | 0% | 100% | 0% |
| **III · animal parasites** (2,395) | 89% | 10% | 85% | 14% | 89% | 10% |
| **IV · animal parasites** (856) | 94% | 6% | 94% | 6% | 94% | 6% |
| IV · plant parasites (4,615) | 93% | 7% | 92% | 7% | 93% | 7% |
| IV · free-living (42) | 100% | 0% | 100% | 0% | 100% | 0% |
| **V · animal parasites** (4,063) | 96% | 4% | 96% | 4% | 96% | 4% |
| V · free-living (326) | 96% | 4% | 96% | 4% | 96% | 4% |
| Other Chromadoria · free-living (15) | 93% | 7% | 93% | 7% | 93% | 7% |
| Family not placed (1) | 0% | 0% | 0% | 0% | 0% | 0% |
<!-- figures:end -->

The figures and tables above are generated from the build by
`bin/readme_figures.py`, alt text included; rerun it after a build to
refresh them in place:

```bash
python3 bin/readme_figures.py work/build docs/figures README.md
```

Groups are the nematode clades of Blaxter et al. (1998), split by
lifestyle. Insect parasites (*Steinernema*, *Heterorhabditis*) count as
animal parasites; fungal feeders (*Aphelenchoides*, *Bursaphelenchus*)
count as plant parasites. A record whose span sits exactly inside
another species' record is tied with it and stops at the rank they share.

Most species-level records are close to another species: 1,019 at one
edit, 836 at two, over a median ITS2 of 240 bases (`species_margins.tsv`).
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
