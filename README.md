# nemabiome-its2-refdb

Builds the nematode ITS2 reference database for NC1/NC2 "nemabiome"
metabarcoding, read on Illumina (2x150, 2x250) and nanopore: every
nematode in the nemabiome ITS2 database, re-cut to the span the reads
cover, graded by how much of it each record holds, with taxonomy,
resolution tables and the adapter look-alike mask. It is the ITS2
counterpart of [vespa-18s-refdb](https://github.com/stephenrdoyle/vespa-18s-refdb),
writes the same release format, and is consumed by
[amphid](https://github.com/stephenrdoyle/nf-amphid-pipeline).

**Status: v0.1.0, the first release (2026-11-10).** See
[docs/DESIGN.md](docs/DESIGN.md) for the design and the measurements
behind it, and [docs/feedback/its2_notes_for_amphid.md](docs/feedback/its2_notes_for_amphid.md)
for what amphid needs to change to use it. DADA2 and IdTaxa users: see
[For DADA2 and IdTaxa users](#for-dada2-and-idtaxa-users).

## Download and verify

```bash
V=v0.1.0
curl -LO https://github.com/stephenrdoyle/nemabiome-its2-refdb/releases/download/$V/nemits2_$V.tar.gz
shasum -a 256 nemits2_$V.tar.gz
# df5a19c2fe043aae09f759eeb2682efe927bdf8ffdbdcfd7fa916cf62ccfe915
tar -xzf nemits2_$V.tar.gz
```

Each file's sha256 is in `nemits2_<version>/manifest.json`; the manifest
and summary of every release are also committed under `releases/`.

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

## v0.1.0

11,956 records from 15,282 source records: 15,108 nemabiome sequences
(version 20260317; 28 dropped as mislabelled and 6 relabelled,
`curation/overrides.tsv`) and 174 ITS2 copies from long-read genome
assemblies and HiFi reads (`curation/long_read_its2.tsv`, which adds six
species: *Loa loa*, *Acanthocheilonema viteae*, *Litomosoides
sigmodontis*, *Spirocerca lupi*, *Trichuris incognita*, *Heligmosomoides
bakeri*). Names follow current NCBI Taxonomy
(`curation/ncbi_names.tsv`; 64 source records renamed). 1,115 partial records are folded into a full-length record of
the same species that holds them. 426 genera, 1,840 species with a
record of their own.

<!-- figures:start -->
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/grade_dark.svg">
  <img alt="Bar chart: records per nematode group in nemits2_v0.1.0, by how much of the ITS2 they hold: I · animal parasites 395 (complete 28, near complete 176, partial 60, unmeasured 131); I · plant parasites 168 (complete 0, near complete 0, partial 0, unmeasured 168); I · free-living 16 (complete 2, near complete 2, partial 0, unmeasured 12); II · plant parasites 56 (complete 0, near complete 0, partial 0, unmeasured 56); II · free-living 19 (complete 1, near complete 0, partial 0, unmeasured 18); III · animal parasites 2,146 (complete 780, near complete 668, partial 484, unmeasured 214); IV · animal parasites 814 (complete 155, near complete 242, partial 111, unmeasured 306); IV · plant parasites 4,299 (complete 2,178, near complete 910, partial 630, unmeasured 581); IV · free-living 42 (complete 13, near complete 9, partial 6, unmeasured 14); V · animal parasites 3,679 (complete 864, near complete 1,609, partial 829, unmeasured 377); V · free-living 306 (complete 143, near complete 62, partial 31, unmeasured 70); Other Chromadoria · free-living 15 (complete 0, near complete 0, partial 0, unmeasured 15); Family not placed 1 (complete 1, near complete 0, partial 0, unmeasured 0)." src="docs/figures/grade_light.svg" width="760">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/taxa_dark.svg">
  <img alt="Bar chart: records per nematode group in nemits2_v0.1.0, and how many resolve to species read whole: I · animal parasites 395 (373 to species); I · plant parasites 168 (168 to species); I · free-living 16 (16 to species); II · plant parasites 56 (56 to species); II · free-living 19 (19 to species); III · animal parasites 2,146 (2,017 to species); IV · animal parasites 814 (775 to species); IV · plant parasites 4,299 (4,080 to species); IV · free-living 42 (42 to species); V · animal parasites 3,679 (3,564 to species); V · free-living 306 (295 to species); Other Chromadoria · free-living 15 (14 to species); Family not placed 1 (0 to species)." src="docs/figures/taxa_light.svg" width="760">
</picture>

How many species nemits2_v0.1.0 can name. A species named in the source records but with no record of its own shares its ITS2 with another species (or sits inside another species' record over its span), so a read from it is named only to genus or above. Names are counted as written.

| Records | Species named in sources | Species with their own record (read whole) | Species named through Illumina 150x150 | Species named through Illumina 250x250 | Genera |
|---|---:|---:|---:|---:|---:|
| All records | 1,885 | 1,840 | 1,833 | 1,840 | 426 |
| **I · animal parasites** | 34 | 33 | 33 | 33 | 6 |
| I · plant parasites | 65 | 65 | 64 | 65 | 3 |
| I · free-living | 10 | 10 | 10 | 10 | 9 |
| II · plant parasites | 18 | 18 | 18 | 18 | 2 |
| II · free-living | 9 | 9 | 9 | 9 | 6 |
| **III · animal parasites** | 299 | 289 | 283 | 289 | 101 |
| **IV · animal parasites** | 152 | 148 | 148 | 148 | 13 |
| IV · plant parasites | 662 | 643 | 643 | 643 | 93 |
| IV · free-living | 21 | 21 | 21 | 21 | 12 |
| **V · animal parasites** | 508 | 503 | 503 | 503 | 156 |
| V · free-living | 104 | 100 | 100 | 100 | 23 |
| Other Chromadoria · free-living | 1 | 1 | 1 | 1 | 2 |
| Family not placed | 2 | 0 | 0 | 0 | 0 |

How far a record can be named in nemits2_v0.1.0, best case (an exact read, no sequencing error): the share of records resolving to species, and to genus at best.

| Records | Nanopore (whole ITS2): species | genus | Illumina 150x150: species | genus | Illumina 250x250: species | genus |
|---|---:|---:|---:|---:|---:|---:|
| All records (11,956) | 96% | 4% | 94% | 6% | 95% | 4% |
| **I · animal parasites** (395) | 94% | 6% | 88% | 12% | 93% | 7% |
| I · plant parasites (168) | 100% | 0% | 99% | 1% | 100% | 0% |
| I · free-living (16) | 100% | 0% | 100% | 0% | 100% | 0% |
| II · plant parasites (56) | 100% | 0% | 100% | 0% | 100% | 0% |
| II · free-living (19) | 100% | 0% | 100% | 0% | 100% | 0% |
| **III · animal parasites** (2,146) | 94% | 6% | 90% | 10% | 94% | 6% |
| **IV · animal parasites** (814) | 95% | 5% | 95% | 5% | 95% | 5% |
| IV · plant parasites (4,299) | 95% | 5% | 94% | 5% | 95% | 5% |
| IV · free-living (42) | 100% | 0% | 100% | 0% | 100% | 0% |
| **V · animal parasites** (3,679) | 97% | 3% | 96% | 4% | 97% | 3% |
| V · free-living (306) | 96% | 2% | 96% | 2% | 96% | 2% |
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

Most species-level records are close to another species: 897 at one
edit, 712 at two, 410 tied (`species_margins.tsv`), over a median ITS2 of
240 bases. At that length amphid's 99.7% species limit allows no
mismatch at all.

## How well reads are named

Tested two ways (docs/DESIGN.md §6): leave-one-out on the database
itself, and two published livestock datasets run through the standard
nemabiome DADA2 workflow.

**Leave-one-out** (3,164 records of 1,210 species, read whole; 2x150 is
within 3 points). *Missing*: every record of the query's species
removed, so any species named is wrong. *Present*: only the query's own
record removed.

| Classifier, this database | Missing: wrong species | Present: right species | Present: wrong species |
|---|---:|---:|---:|
| `bin/classify_asvs.py`, gap rule | **6.5%** | 52.6% | **0.7%** |
| DECIPHER IdTaxa (threshold 60) | 22.0% | 77.1% | 2.5% |
| dada2 assignTaxonomy (minBoot 80) | 42.4% | **85.1%** | 4.6% |

The gap rule names a species only when the read is >= 99% identical to
it and every other species is >= 4 edits further off; otherwise it stops
at genus (>= 97%) or family (>= 95%). It is the most conservative of the
three: it names fewest species but is wrong least, which matters most in
samples with species missing from the database (wildlife, new hosts).
Tuning the rule per genus or family did not help on held-out species
(`results/leave_one_out/group_rules.txt`). Each call carries `nearest`
(the best species and its identity) and `next` (the nearest other
species and how many edits further), so a genus call still says how
close it came: *Trichostrongylus*; nearest *T. retortaeformis* 100.0%;
next *T. axei* +2.

**Real data** (share of reads named to species; all arms named the same
species wherever two of them named one):

| Dataset | assignTaxonomy x nemabiome | IdTaxa x nemabiome | assignTaxonomy x this db | IdTaxa x this db | gap rule x this db |
|---|---:|---:|---:|---:|---:|
| Šlapeta sheep, PRJNA936713, 2x250 | 69.5% | 95.4% | 99.3% | 99.7% | 84.1% |
| same, cut to 2x150 | 69.1% | 87.9% | 99.2% | 92.2% | 52.4% |
| Gilleard cattle, PRJNA1524548, 2x250 | 72.6% | 68.0% | 92.9% | 69.9% | 62.4% |
| same, cut to 2x150 | 69.0% | 67.8% | 93.7% | 91.9% | 37.3% |

assignTaxonomy and IdTaxa compare k-mers over the whole ASV, so with
this database the ASVs' 28S stretch must be cut off first
(`bin/trim_28s.py`); the columns for this database are on cut ASVs.
The gap rule's figures are lower than those first published with v0.1.0
(92.2 / 90.5 / 86.9 / 86.9%): the classifier then scored only the top 40
candidates and could miss the nearest other species (issue #1). Most of
what it now stops at genus is *Haemonchus*, held back by a few partial
single-accession records 2-3 edits from *H. placei* and *H. contortus*
(issue #2; docs/DESIGN.md §6.2).
`results/dada2_comparison/` holds each arm's summary and the per-ASV
calls.

**Simulated mocks** (Jesudoss Chelladurai et al. 2025; canine, equine,
ruminant, all species present): no classifier, with either database,
named a wrong species for an exact copy of a source sequence; the
k-mer classifiers named 86-100% of reads to species, the gap rule
45-100%, stopping at genus where a sister species is 1-3 edits off.
The mocks' truth includes two backbone mislabels this database corrects
(docs/DESIGN.md §6.3; `results/mocks/`).

## For DADA2 and IdTaxa users

The release carries the database in the formats of the nemabiome.ca
DADA2 workflow, to use in place of nemabiome's own files:

| file | use |
|---|---|
| `nemits2_<version>_dada2.fasta` | `assignTaxonomy` reference |
| `nemits2_<version>_idtaxa.RData` | `load()` it for `trainingSet`, then `IdTaxa` |
| `trim_28s.R` | **required**: cut the 28S stretch off ASVs before classifying |

Records here stop where the ~47-base 28S stretch before NC2 starts; your
ASVs end in it. assignTaxonomy and IdTaxa compare k-mers over the whole
ASV, so uncut ASVs lose most species calls (one sheep set: 45% of reads
to species uncut, 99.3% cut). After `removeBimeraDenovo`:

```r
source("trim_28s.R")
seqs <- getSequences(seqtab.nochim)
cut  <- trim_28s(seqs)
ok   <- nchar(gsub("N", "", cut)) >= 100      # shorter is not ITS2

taxa <- assignTaxonomy(cut[ok], "nemits2_<version>_dada2.fasta", minBoot = 80)
rownames(taxa) <- seqs[ok]                      # back to your ASVs

load("nemits2_<version>_idtaxa.RData")          # trainingSet
ids <- IdTaxa(DNAStringSet(cut[ok]), trainingSet, strand = "both", threshold = 60)
```

Concatenated 2x150 ASVs (`mergePairs(justConcatenate = TRUE)`) are cut
in their reverse part. Set a seed and report it: IdTaxa's calls are
random draws, and ASVs near its threshold change between runs - in one
cattle set, *Haemonchus placei* ASVs holding 10-23% of reads sat at
confidence 62-76 and switched between species and genus with the seed.

Which classifier to use is a trade-off (docs/DESIGN.md §6.2): k-mer
classifiers name the most species when every species is in the
database; when some may be missing, they name a wrong one far more often
than the span-scored gap rule (`bin/classify_asvs.py`).

## Build

```bash
S=sources/nemabiome_its2_20260317/Nematode_ITS2_20260317
python3 bin/nemabiome_export.py sources/nemabiome_its2_20260317/raw/database.xlsx $S
python3 bin/fetch_backbone_records.py ${S}_taxonomy_nr.txt sources/genbank_backbone/records.fasta.gz
python3 bin/its2.py ${S}_seqs_nr.fasta sources/genbank_backbone/records.fasta.gz work/its2/its2_cut
python3 bin/its2_partials.py work/its2/its2_cut work/its2/its2_refs 8
python3 bin/ncbi_names.py ${S}_taxonomy_nr.txt curation/ncbi_names.tsv   # current NCBI names; refetch before a release
python3 bin/build.py work/its2/its2 ${S}_taxonomy_nr.txt nemabiome_20260317 work/build --workers 8 \
    --layer curation/long_read_its2.tsv
bin/export_kmer.sh work/build work/build/kmer nemits2_dev   # DADA2 / IdTaxa files (R, dada2, DECIPHER)
```

To cut a release, from committed code: build with `--version vX.Y.Z` into
a fresh directory, then `bin/release.py <build dir> releases vX.Y.Z`. It
adds the DADA2 / IdTaxa files, re-checksums, and writes
`releases/nemits2_vX.Y.Z/` and its tarball; it refuses a build from
uncommitted code and never overwrites a release.

The long-read layer is committed. To remake it, scan each assembly or
read set listed in `curation/long_read_sources.tsv` with
`bin/its2_scan.py scan`, then
`bin/genome_layer.py curation/long_read_sources.tsv ${S}_taxonomy_nr.txt curation/long_read_its2.tsv`.

Testing: `bin/its2_leave_one_out.py` (gap rule), `bin/loo_dada2.py`
(assignTaxonomy or IdTaxa on the same queries), `bin/loo_group_rules.py`;
`bin/dada2_nemabiome.R`, `bin/trim_28s.py`, `bin/assign_taxonomy.R`,
`bin/idtaxa.R`, `bin/classify_asvs.py` and `bin/compare_arms.py` for
real data; `bin/mock_arms.sh` and `bin/score_mocks.py` for the mocks. These need R with dada2 and DECIPHER, and cutadapt.

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
| `nemits2_<version>_dada2.fasta`, `nemits2_<version>_idtaxa.RData`, `trim_28s.R` | for the DADA2 workflow (above) |

## Citation

The backbone is the nemabiome ITS2 database (Workentine et al. 2020,
BMC Genetics 21:74; https://www.nemabiome.ca), CC BY 4.0. Primers:
Gasser et al. 1993, Nucleic Acids Res 21:2525; method: Avramenko et al.
2015, PLoS ONE 10:e0143559.

MIT licence for the code (LICENSE).
