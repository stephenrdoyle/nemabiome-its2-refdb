# nemabiome-its2-refdb

Builds the nematode ITS2 reference database for NC1/NC2 "nemabiome"
metabarcoding, read on Illumina (2x150, 2x250) and nanopore: every
nematode in the nemabiome ITS2 database, re-cut to the span the reads
cover, graded by how much of it each record holds, with taxonomy,
resolution tables and the adapter look-alike mask. It is the ITS2
counterpart of [vespa-18s-refdb](https://github.com/stephenrdoyle/vespa-18s-refdb),
writes the same release format, and is consumed by
[amphid](https://github.com/stephenrdoyle/nf-amphid-pipeline).

**Status: v0.1.3 (2026-11-11).** v0.1.0 was the first release. See
[docs/DESIGN.md](docs/DESIGN.md) for the design and the measurements
behind it, and [docs/feedback/its2_notes_for_amphid.md](docs/feedback/its2_notes_for_amphid.md)
for what amphid needs to change to use it. DADA2 and IdTaxa users: see
[For DADA2 and IdTaxa users](#for-dada2-and-idtaxa-users).

## Download and verify

```bash
V=v0.1.3
curl -LO https://github.com/stephenrdoyle/nemabiome-its2-refdb/releases/download/$V/nemits2_$V.tar.gz
shasum -a 256 nemits2_$V.tar.gz
# 890beb48998fa6a1f2203e4184b342a380f2274de663a73dd5a4a1f6c9b9385e
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

## v0.1.3

v0.1.3 adds a `lineage` column to `taxonomy.tsv` (the *Teladorsagia
circumcincta* goat line, `curation/lineages.tsv`); records and calls are
as in v0.1.2.

11,135 records from 15,282 source records: 15,108 nemabiome sequences
(version 20260317; 54 dropped and 33 relabelled,
`curation/overrides.tsv`) and 174 ITS2 copies from long-read genome
assemblies and HiFi reads (`curation/long_read_its2.tsv`, which adds six
species: *Loa loa*, *Acanthocheilonema viteae*, *Litomosoides
sigmodontis*, *Spirocerca lupi*, *Trichuris incognita*, *Heligmosomoides
bakeri*). Names follow current NCBI Taxonomy
(`curation/ncbi_names.tsv`; 64 source records renamed). 1,125 partial
records are folded into a full-length record of the same species that
holds them, and 749 weak records (partial or unmeasured, one accession)
into a record of their species within 1 edit. 426 genera, 1,836 species
with a record of their own.

Changes from v0.1.1 (11,146 records): the 2026-11-11 curation pass
(`docs/curation_pass_2026-11-11.md`). 27 relabels fold minor male morphs
into their species (*Cooperia surnabada* and *C. spatulata*,
*Marshallagia occidentalis*, *Teladorsagia davtiani*, *Ostertagia
arctica*); 26 drops remove records that sit next to another species
while the rest of their own sits far off, and shared ITS2s typical of
only one of their names. Leave-one-out species recall rises from 53.4%
to 55.3% and over-classification falls from 6.4% to 5.9%; *Cooperia*
reads are now named to species. v0.1.1 changed from v0.1.0 (11,956
records): weak records folded (issue #2), 36 records deposited
reverse-complemented turned and re-cut, primer and 28S remnants trimmed,
and margins and resolution tables computed against every candidate
(issue #1).

<!-- figures:start -->
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/grade_dark.svg">
  <img alt="Bar chart: records per nematode group in nemits2_v0.1.3, by how much of the ITS2 they hold: I · animal parasites 390 (complete 28, near complete 176, partial 55, unmeasured 131); I · plant parasites 168 (complete 0, near complete 0, partial 0, unmeasured 168); I · free-living 16 (complete 2, near complete 2, partial 0, unmeasured 12); II · plant parasites 56 (complete 0, near complete 0, partial 0, unmeasured 56); II · free-living 19 (complete 1, near complete 0, partial 0, unmeasured 18); III · animal parasites 2,003 (complete 778, near complete 670, partial 348, unmeasured 207); IV · animal parasites 702 (complete 155, near complete 240, partial 89, unmeasured 218); IV · plant parasites 4,098 (complete 2,186, near complete 899, partial 450, unmeasured 563); IV · free-living 42 (complete 13, near complete 9, partial 6, unmeasured 14); V · animal parasites 3,334 (complete 1,077, near complete 1,431, partial 557, unmeasured 269); V · free-living 291 (complete 173, near complete 58, partial 19, unmeasured 41); Other Chromadoria · free-living 15 (complete 0, near complete 0, partial 0, unmeasured 15); Family not placed 1 (complete 1, near complete 0, partial 0, unmeasured 0)." src="docs/figures/grade_light.svg" width="760">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/taxa_dark.svg">
  <img alt="Bar chart: records per nematode group in nemits2_v0.1.3, and how many resolve to species read whole: I · animal parasites 390 (368 to species); I · plant parasites 168 (168 to species); I · free-living 16 (16 to species); II · plant parasites 56 (56 to species); II · free-living 19 (19 to species); III · animal parasites 2,003 (1,880 to species); IV · animal parasites 702 (684 to species); IV · plant parasites 4,098 (3,938 to species); IV · free-living 42 (42 to species); V · animal parasites 3,334 (3,261 to species); V · free-living 291 (280 to species); Other Chromadoria · free-living 15 (14 to species); Family not placed 1 (0 to species)." src="docs/figures/taxa_light.svg" width="760">
</picture>

How many species nemits2_v0.1.3 can name. A species named in the source records but with no record of its own shares its ITS2 with another species (or sits inside another species' record over its span), so a read from it is named only to genus or above. Names are counted as written.

| Records | Species named in sources | Species with their own record (read whole) | Species named through Illumina 150x150 | Species named through Illumina 250x250 | Genera |
|---|---:|---:|---:|---:|---:|
| All records | 1,880 | 1,836 | 1,829 | 1,836 | 426 |
| **I · animal parasites** | 34 | 33 | 33 | 33 | 6 |
| I · plant parasites | 65 | 65 | 64 | 65 | 3 |
| I · free-living | 10 | 10 | 10 | 10 | 9 |
| II · plant parasites | 18 | 18 | 18 | 18 | 2 |
| II · free-living | 9 | 9 | 9 | 9 | 6 |
| **III · animal parasites** | 299 | 289 | 283 | 289 | 101 |
| **IV · animal parasites** | 152 | 148 | 148 | 148 | 13 |
| IV · plant parasites | 662 | 643 | 643 | 643 | 93 |
| IV · free-living | 21 | 21 | 21 | 21 | 12 |
| **V · animal parasites** | 503 | 499 | 499 | 499 | 156 |
| V · free-living | 104 | 100 | 100 | 100 | 23 |
| Other Chromadoria · free-living | 1 | 1 | 1 | 1 | 2 |
| Family not placed | 2 | 0 | 0 | 0 | 0 |

How far a record can be named in nemits2_v0.1.3, best case (an exact read, no sequencing error): the share of records resolving to species, and to genus at best.

| Records | Nanopore (whole ITS2): species | genus | Illumina 150x150: species | genus | Illumina 250x250: species | genus |
|---|---:|---:|---:|---:|---:|---:|
| All records (11,135) | 96% | 3% | 95% | 5% | 96% | 4% |
| **I · animal parasites** (390) | 94% | 6% | 88% | 12% | 93% | 7% |
| I · plant parasites (168) | 100% | 0% | 99% | 1% | 100% | 0% |
| I · free-living (16) | 100% | 0% | 100% | 0% | 100% | 0% |
| II · plant parasites (56) | 100% | 0% | 100% | 0% | 100% | 0% |
| II · free-living (19) | 100% | 0% | 100% | 0% | 100% | 0% |
| **III · animal parasites** (2,003) | 94% | 6% | 90% | 10% | 94% | 6% |
| **IV · animal parasites** (702) | 97% | 3% | 97% | 3% | 97% | 3% |
| IV · plant parasites (4,098) | 96% | 4% | 96% | 4% | 96% | 4% |
| IV · free-living (42) | 100% | 0% | 100% | 0% | 100% | 0% |
| **V · animal parasites** (3,334) | 98% | 2% | 97% | 3% | 98% | 2% |
| V · free-living (291) | 96% | 2% | 96% | 2% | 96% | 2% |
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
itself, and two livestock datasets (one published) run through the standard
nemabiome DADA2 workflow.

Figures in this section are for v0.1.2 and v0.1.3 (same records). Since v0.1.0 the classifier
also scores every candidate (issue #1) and scores read windows strictly
(docs/DESIGN.md §6).

**Leave-one-out** (3,182 records of 1,220 species, read whole; 2x150 is
within 3 points). *Present*: only the query's own record removed.
*Missing*: every record of the query's species removed, so any species
named is wrong. Recall = right calls / reads; precision = right calls /
calls made at that rank; over-classification = species named for a
missing species (`results/leave_one_out/metrics.txt`).

| Classifier, this database | Species recall | Species precision | Genus recall | Genus precision | Over-classification |
|---|---:|---:|---:|---:|---:|
| `bin/classify_asvs.py`, gap rule | 55.3% | **98.8%** | 88.9% | 99.3% | **5.9%** |
| DECIPHER IdTaxa (threshold 60) | 78.8% | 96.9% | 93.5% | **99.5%** | 23.2% |
| dada2 assignTaxonomy (minBoot 80) | **86.0%** | 95.1% | **97.2%** | 99.2% | 43.5% |

The gap rule names a species only when the read is >= 99% identical to
it and every other species is >= 4 edits further off, or the match is
exact and every other species is >= 2 edits off; otherwise it stops at
genus (>= 97%) or family (>= 95%). It is the most conservative of the
three: it names fewest species but is wrong least, which matters most in
samples with species missing from the database (wildlife, new hosts). At
genus all three are about 99% precise; they differ in how often they
name a species that is not there.
Tuning the rule per genus or family did not help on held-out species
(`results/leave_one_out/group_rules.txt`). Each call carries `nearest`
(the best species and its identity) and `next` (the nearest other
species and how many edits further), so a genus call still says how
close it came: *Trichostrongylus*; nearest *T. colubriformis* 100.0%;
next *T. capricola* +1.

**Real data** (share of reads named to species; all arms named the same
species wherever two of them named one). Two NC1/NC2 sets from the
Šlapeta group, one published,
2x250 pairs joined hybrid (merged where they overlap, concatenated where
they do not), and cut to 2x150:

| Dataset | assignTaxonomy x nemabiome | IdTaxa x nemabiome | assignTaxonomy x this db | IdTaxa x this db | gap rule x this db |
|---|---:|---:|---:|---:|---:|
| Francis & Šlapeta, cattle, goats, sheep, PRJNA936713, 2x250 | 72.0% | 81.7% | 99.9% | 99.8% | 93.8% |
| same, cut to 2x150 | 76.0% | 88.2% | 99.1% | 99.1% | 94.5% |
| Šlapeta group, unpublished, PRJNA1524548, cattle, kangaroos, pigs, rabbits, 2x250 | 70.4% | 69.7% | 94.9% | 99.1% | 95.2% |
| same, cut to 2x150 | 73.7% | 65.9% | 94.8% | 87.2% | 95.9% |

assignTaxonomy and IdTaxa compare k-mers over the whole ASV, so with
this database the ASVs' 28S stretch must be cut off first
(`bin/trim_28s.py`); the columns for this database are on cut ASVs.
Merging only, as the standard workflow does, loses the ITS2 too long
for a 2x250 pair to overlap: two of the four PRJNA1524548 pig libraries kept
6 and 665 reads merged and 22,932 and 19,801 hybrid, all *Metastrongylus
apri* (see For DADA2 and IdTaxa users). Half the PRJNA1524548 libraries have
read 2 starting inside 28S, without NC2; `bin/trim_primers.py` keeps
them. `results/dada2_comparison/` holds each arm's summary, the per-ASV
calls, the merged and concatenated runs, and an HTML report per set
(docs/DESIGN.md §6.2).

**Simulated mocks** (Jesudoss Chelladurai et al. 2025; canine, equine,
ruminant, all species present, so they measure recall, not
over-classification): no classifier, with either database, named a
wrong species for an exact copy of a source sequence (precision 100%);
the k-mer classifiers named 89-100% of reads to species, the gap rule
88-100%, stopping at genus where a sister species is 1 edit off.
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
ASV, so uncut ASVs lose most species calls (one livestock set: 45% of
reads to species uncut, 99.3% cut).

At 2x250, join pairs hybrid rather than by merging alone. `mergePairs`
needs 12 bases of overlap; after primer trimming (about 230 + 230
bases) an ITS2 over about 400 bases never merges, and those pairs are
silently dropped: *Metastrongylus* (413-422 bases), Trichuridae,
Metastrongylidae, Onchocercidae and others. Keep the merged pairs and
concatenate only the rest:

```r
mm <- mergePairs(dadaF, derepF, dadaR, derepR, returnRejects = TRUE)
cc <- mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate = TRUE)
hybrid <- mapply(function(m, c) {
  acc <- m[m$accept, ]
  rest <- !paste(c$forward, c$reverse) %in% paste(acc$forward, acc$reverse)
  rbind(acc[, c("sequence", "abundance")], c[rest, c("sequence", "abundance")])
}, mm, cc, SIMPLIFY = FALSE)
seqtab <- makeSequenceTable(hybrid)
```

(for one sample, `mm` and `cc` are data frames rather than lists: wrap
them in `list()`; `bin/dada2_nemabiome.R` does this as its `hybrid`
mode). Concatenating every pair instead suits the k-mer classifiers but
costs the gap rule species calls where read 1 runs through a short ITS2
into 28S. After `removeBimeraDenovo`:

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

Concatenated ASVs (2x150, or the unmerged part of a hybrid join) are
cut in their reverse part. Set a seed and report it: IdTaxa's calls are
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
`bin/real_data.sh` for real data (`bin/trim_primers.py`,
`bin/dada2_nemabiome.R`, `bin/trim_28s.py`, `bin/assign_taxonomy.R`,
`bin/idtaxa.R`, `bin/classify_asvs.py`, `bin/compare_arms.py`, then
`bin/real_data_report.py` for the HTML report); `bin/mock_arms.sh` and
`bin/score_mocks.py` for the mocks. These need R with dada2 and DECIPHER.

Needs Python 3.9+ and minimap2. `sources/README.md` says where each
download comes from.

## Release files

| file | holds |
|---|---|
| `refdb.fasta` | one record per unique ITS2 |
| `taxonomy.tsv` | superkingdom ... species, deepest rank, members, grade, coverage, sources, lineage (a published line within a species; blank for nearly all) |
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
