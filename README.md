# nematode-its2-refdb

Builds a curated nematode ITS2 reference database for NC1/NC2
"nemabiome" metabarcoding, read on Illumina (2x150, 2x250) and
nanopore. It starts from the nemabiome ITS2 database and adds GenBank
records published since, informal "Genus sp." lineages, and ITS2 copies
from genome assemblies; each record is re-cut to the span the reads
cover, graded by how much of it it holds, and screened, with mislabelled
records dropped or relabelled (`curation/`). It writes taxonomy,
resolution tables and the adapter look-alike mask. It is an independent
derivative, not a product of the nemabiome project; see
[Citation](#citation). It is the ITS2
counterpart of [vespa-18s-refdb](https://github.com/stephenrdoyle/vespa-18s-refdb),
writes the same release format, and is consumed by
[amphid](https://github.com/stephenrdoyle/nf-amphid-pipeline).

**Status: v0.1.4 (2026-11-12).** v0.1.0 was the first release. See
[docs/DESIGN.md](docs/DESIGN.md) for the design and the measurements
behind it, and [docs/feedback/its2_notes_for_amphid.md](docs/feedback/its2_notes_for_amphid.md)
for what amphid needs to change to use it. DADA2 and IdTaxa users: see
[For DADA2 and IdTaxa users](#for-dada2-and-idtaxa-users).

## Download and verify

```bash
V=v0.1.4
curl -LO https://github.com/stephenrdoyle/nematode-its2-refdb/releases/download/$V/nemits2_$V.tar.gz
shasum -a 256 nemits2_$V.tar.gz
# e8ed2dc74a5010a5a6f621329ad6f0a34da3cf9c794c3091a9d3dc2f78f165d6
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

## Since v0.1.4 (on `main`, not yet released)

Sixteen more genome assemblies scanned for ITS2 (`curation/long_read_sources.tsv`;
docs/DESIGN.md §5.5): 74 copies of 15 species, 49 of them new records
(*Dirofilaria immitis*, *Acanthocheilonema viteae*, *Dictyocaulus
viviparus*, *Trichuris trichiura*, *Thelazia callipaeda* and others).
The genomes settle two curation questions (`curation/overrides.tsv`):

- **The only *Trichostrongylus capricola* record, JF276023, is dropped.**
  It is 1 edit from 692 of the 766 ITS2 copies of the
  *T. colubriformis* genome GCA_963978875, inside that species' own
  spread. *T. colubriformis* reads it held at genus are now named:
  5.8% of the ruminant mock's reads (gap rule 90.8% -> 96.9% to
  species), 2,052 Francis & Šlapeta reads (goats 87.2% -> 88.6%).
- **Eight Sri Lankan "*D. repens*" records are *D. asiatica***
  (ON332065-71, PV643877): 1-3 edits from the *D. asiatica* genome,
  14-21 from the *D. repens* genome.

Leave-one-out (3,230 queries): present species named 57.3% (gap rule) /
68.2% (tiered), as in v0.1.4; over-classification 5.4% / 9.0%, from
5.5% / 9.3%. The mocks still name no wrong species
(`results/`, docs/DESIGN.md §6).

***Metastrongylus* curation** (docs/DESIGN.md §8.18;
`curation/genbank_named.tsv`, `curation/overrides.tsv`). Most records
are SSCP bands of Conole et al. (2001, Estonian wild boar), whose ITS2
varies chiefly in microsatellites. *M. pudendotectus* is one species
(records up to 78 edits apart plain are 2-11 apart with repeats
collapsed, and repeat-aware scoring names it); *M. salmi* is not: three
tight types, A and C 21-34 edits apart, about as far as A is from
*M. apri*. Type C is the sequence Humbert & Henry (1997) deposited as
*M. asymmetricus*, so it stays tied with that species at genus. Type B
(AJ305395/400/403) is type A lacking a 44-base block every other record
carries, most likely a gel-band artefact, and is dropped. Added by hand:
the rest of the Conole bands, *M. apri* from Denmark, Viet Nam and
Brazil, and 14 Brazilian wild boar *M. salmi* (OR660330-50, unpublished),
a further type 17-27 edits from every other record.

**Informal "Genus sp." records** (docs/DESIGN.md §5.7). Some lineages
are known only by an informal name, such as the gorilla hookworm
"*Necator* sp." of Hasegawa et al. The GenBank update now takes every
nematode ITS2 record named "Genus sp. X"; the build gives each its genus
only (`species_status` `unnamed`, the name in `lineage`) and drops it
where a named record's ITS2 is identical to, inside or contains it, or
where its genus disagrees with its ITS2 (a named record of another
genus >= 97% identical and closer than any of its own genus: 86 records,
mostly *Raphidascaris* and *Terranova* inside *Hysterothylacium* and
*Euterranova*). The classifier scores reads against named records
first and uses informal ones only where that call stops above genus, so
they add genus calls and never change one. 737 are kept. Leave-one-out:
missing species reach the right genus 23.5% of the time (from 19.3%),
wrong genus unchanged at 1.2%. Real data: no species call changes;
gorilla *Necator* (18-21% of Cameroon ape reads) and kangaroo
*Macropostrongylus* (0.6-0.7% of PRJNA1524548) reach genus instead of
unassigned.

The development build has 12,054 records and 1,873 named species.

## v0.1.4

11,288 records and 1,874 named species (427 genera) from 16,051 source
records: the nemabiome backbone plus 733 GenBank records it lacks
(`bin/genbank_update.py`: ITS2 records published since the backbone's
17 March 2026 snapshot, 5.8S/28S records holding an ITS2, and 98
accessions listed by hand in `curation/genbank_named.tsv`). Update
records are taken only where `bin/its2.py` cut them at an anchor or
their own end; 1,057 were not. New species include *Heligmosomoides
bibullosus*, *Metastrongylus asymmetricus* and *M. confusus*.

- **Informal names are not species.** "sp."/"spp." names (*Oesophagostomum*
  sp. 3335_MDR_CAM_OH15, *Cloacina* sp. 1 AJ-2019, 19 records) are
  genus-level, unnamed records; the informal name is in the `lineage`
  column (`curation/lineages.tsv`).
- **Species coherence screen** (`bin/species_coherence.py`,
  `results/curation/species_coherence.tsv`): records far from their own
  species and next to another. 65 more drops and 2 relabels
  (`curation/overrides.tsv`, `docs/curation_pass_2026-11-11.md` F-H),
  among them *Ashworthius sidemi* OP787564 (an *H. contortus* ITS2) and
  24 West African "*H. contortus*" carrying the *H. placei* ITS2.
- **Cuts:** an NC1 well inside a cut is taken as the start (5.8S and
  ITS1 no longer ride along), and a record ending in a short 28S stub is
  cut at it (a complete *Oesophagostomum radiatum*, PZ347037, is now in).
- **Repeat-aware scoring** (`bin/its2_classify.py`): a read that differs
  from a reference only in how many units a simple repeat has is one edit
  off, not one per base. Goat *Oe. venulosum* reads with a long (AT)n
  stopped at family; they are now named (docs/DESIGN.md §6).
- **A tentative tier** (`tiered_rank`, `tiered_name`, `tentative` in
  `bin/classify_asvs.py`'s output): species at >= 98% with the next
  species >= 4 edits off, flagged, beside the gap rule's confident call.
- `bin/dada2_nemabiome.R` takes `TRUNCLEN="F,R"` for a fixed truncation
  (tested; it cost species calls, so the default stays off).

### Earlier releases

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
  <img alt="Bar chart: records per nematode group in nemits2_v0.1.4, by how much of the ITS2 they hold: I · animal parasites 371 (complete 28, near complete 176, partial 55, unmeasured 112); I · plant parasites 167 (complete 0, near complete 0, partial 0, unmeasured 167); I · free-living 16 (complete 2, near complete 2, partial 0, unmeasured 12); II · plant parasites 56 (complete 0, near complete 0, partial 0, unmeasured 56); II · free-living 19 (complete 1, near complete 0, partial 0, unmeasured 18); III · animal parasites 2,016 (complete 786, near complete 667, partial 356, unmeasured 207); IV · animal parasites 701 (complete 156, near complete 242, partial 91, unmeasured 212); IV · plant parasites 4,165 (complete 2,205, near complete 906, partial 480, unmeasured 574); IV · free-living 42 (complete 13, near complete 9, partial 6, unmeasured 14); V · animal parasites 3,416 (complete 1,113, near complete 1,453, partial 589, unmeasured 261); V · free-living 303 (complete 199, near complete 47, partial 20, unmeasured 37); Other Chromadoria · free-living 15 (complete 0, near complete 0, partial 0, unmeasured 15); Family not placed 1 (complete 1, near complete 0, partial 0, unmeasured 0)." src="docs/figures/grade_light.svg" width="760">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/figures/taxa_dark.svg">
  <img alt="Bar chart: records per nematode group in nemits2_v0.1.4, and how many resolve to species read whole: I · animal parasites 371 (365 to species); I · plant parasites 167 (167 to species); I · free-living 16 (16 to species); II · plant parasites 56 (56 to species); II · free-living 19 (19 to species); III · animal parasites 2,016 (1,892 to species); IV · animal parasites 701 (691 to species); IV · plant parasites 4,165 (3,982 to species); IV · free-living 42 (42 to species); V · animal parasites 3,416 (3,324 to species); V · free-living 303 (292 to species); Other Chromadoria · free-living 15 (14 to species); Family not placed 1 (0 to species)." src="docs/figures/taxa_light.svg" width="760">
</picture>

How many species nemits2_v0.1.4 can name. A species named in the source records but with no record of its own shares its ITS2 with another species (or sits inside another species' record over its span), so a read from it is named only to genus or above. Names are counted as written.

| Records | Species named in sources | Species with their own record (read whole) | Species named through Illumina 150x150 | Species named through Illumina 250x250 | Genera |
|---|---:|---:|---:|---:|---:|
| All records | 1,904 | 1,846 | 1,840 | 1,846 | 427 |
| **I · animal parasites** | 34 | 33 | 33 | 33 | 6 |
| I · plant parasites | 65 | 65 | 64 | 65 | 3 |
| I · free-living | 10 | 10 | 10 | 10 | 9 |
| II · plant parasites | 18 | 18 | 18 | 18 | 2 |
| II · free-living | 9 | 9 | 9 | 9 | 6 |
| **III · animal parasites** | 300 | 290 | 285 | 290 | 101 |
| **IV · animal parasites** | 152 | 149 | 149 | 149 | 13 |
| IV · plant parasites | 665 | 646 | 646 | 646 | 94 |
| IV · free-living | 21 | 21 | 21 | 21 | 12 |
| **V · animal parasites** | 521 | 502 | 502 | 502 | 156 |
| V · free-living | 106 | 102 | 102 | 102 | 23 |
| Other Chromadoria · free-living | 1 | 1 | 1 | 1 | 2 |
| Family not placed | 2 | 0 | 0 | 0 | 0 |

How far a record can be named in nemits2_v0.1.4, best case (an exact read, no sequencing error): the share of records resolving to species, and to genus at best.

| Records | Nanopore (whole ITS2): species | genus | Illumina 150x150: species | genus | Illumina 250x250: species | genus |
|---|---:|---:|---:|---:|---:|---:|
| All records (11,288) | 96% | 4% | 95% | 5% | 96% | 4% |
| **I · animal parasites** (371) | 98% | 2% | 92% | 8% | 97% | 3% |
| I · plant parasites (167) | 100% | 0% | 99% | 1% | 100% | 0% |
| I · free-living (16) | 100% | 0% | 100% | 0% | 100% | 0% |
| II · plant parasites (56) | 100% | 0% | 100% | 0% | 100% | 0% |
| II · free-living (19) | 100% | 0% | 100% | 0% | 100% | 0% |
| **III · animal parasites** (2,016) | 94% | 6% | 90% | 10% | 94% | 6% |
| **IV · animal parasites** (701) | 99% | 1% | 99% | 1% | 99% | 1% |
| IV · plant parasites (4,165) | 96% | 4% | 95% | 5% | 96% | 4% |
| IV · free-living (42) | 100% | 0% | 100% | 0% | 100% | 0% |
| **V · animal parasites** (3,416) | 97% | 3% | 97% | 3% | 97% | 3% |
| V · free-living (303) | 96% | 2% | 96% | 2% | 96% | 2% |
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

Figures in this section are for v0.1.4. Since v0.1.0 the classifier
also scores every candidate (issue #1) and scores read windows strictly
(docs/DESIGN.md §6).

**Leave-one-out** (3,224 records of 1,236 species, read whole; 2x150 is
within 3 points). *Present*: only the query's own record removed.
*Missing*: every record of the query's species removed, so any species
named is wrong. Recall = right calls / reads; precision = right calls /
calls made at that rank; over-classification = species named for a
missing species (`results/leave_one_out/metrics.txt`).

| Classifier, this database | Species recall | Species precision | Genus recall | Genus precision | Over-classification |
|---|---:|---:|---:|---:|---:|
| `bin/classify_asvs.py`, gap rule | 57.3% | **99.2%** | 89.6% | 99.4% | **5.5%** |
| `bin/classify_asvs.py`, tiered (gap + tentative) | 68.2% | 98.9% | 89.6% | 99.4% | 9.3% |
| DECIPHER IdTaxa (threshold 60) | 80.3% | 97.5% | 94.5% | **99.6%** | 20.6% |
| dada2 assignTaxonomy (minBoot 80) | **86.3%** | 95.4% | **97.7%** | 99.3% | 42.9% |

The gap rule names a species only when the read is >= 99% identical to
it and every other species is >= 4 edits further off, or the match is
exact and every other species is >= 2 edits off; otherwise it stops at
genus (>= 97%) or family (>= 95%). It is the most conservative of the
three: it names fewest species but is wrong least, which matters most in
samples with species missing from the database (wildlife, new hosts). At
genus all three are about 99% precise; they differ in how often they
name a species that is not there. Differences in a simple repeat's
length count as one edit when finding the best match, but not when
measuring the gap to the next species.

The tiered rule also names, as *tentative*, a species the read is
>= 98% identical to with every other species >= 4 edits off: reads
1-2% from their own species (the *Teladorsagia circumcincta* goat line,
*Oe. venulosum*), but also reads of a missing sister species. Read
tentative calls as species where the hosts' worms are well sampled
(livestock), as genus where they may not be (wildlife).
Tuning the rule per genus or family did not help on held-out species
(`results/leave_one_out/group_rules.txt`). Each call carries `nearest`
(the best species and its identity) and `next` (the nearest other
species and how many edits further), so a genus call still says how
close it came: *Trichostrongylus*; nearest *T. rugatus* 99.6%;
next *T. retortaeformis* +2.

**Real data** (share of reads named to species; all arms named the same
species wherever two of them named one). Two NC1/NC2 sets from the
Šlapeta group, one published,
2x250 pairs joined hybrid (merged where they overlap, concatenated where
they do not), and cut to 2x150:

| Dataset | assignTaxonomy x nemabiome | IdTaxa x nemabiome | assignTaxonomy x this db | IdTaxa x this db | gap rule x this db | tiered x this db |
|---|---:|---:|---:|---:|---:|---:|
| Francis & Šlapeta, cattle, goats, sheep, PRJNA936713, 2x250 | 72.0% | 81.7% | 99.5% | 99.8% | 95.3% | 99.0% |
| same, cut to 2x150 | 76.0% | 88.2% | 99.1% | 99.7% | 94.8% | 98.1% |
| Šlapeta group, unpublished, PRJNA1524548, cattle, kangaroos, pigs, rabbits, 2x250 | 70.4% | 69.7% | 95.1% | 75.9%* | 94.3% | 98.3% |
| same, cut to 2x150 | 73.7% | 65.9% | 94.8% | 85.5%* | 94.7% | 98.3% |

\* IdTaxa with seed 1, which leaves the main *H. placei* ASVs at genus;
seeds 2 and 3 name 99.0% and 99.2% of the 2x250 reads to species.

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
the k-mer classifiers named 88-100% of reads to species, the gap rule
91-100%, stopping at genus where a sister species is 1 edit off; the
tiered rule named no wrong species either.
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
than the span-scored gap rule (`bin/classify_asvs.py`). Its tiered rule
sits between them: 98-99% of reads named on both livestock sets, with
less than half IdTaxa's over-classification.

## Build

```bash
S=sources/nemabiome_its2_20260317/Nematode_ITS2_20260317
python3 bin/nemabiome_export.py sources/nemabiome_its2_20260317/raw/database.xlsx $S
python3 bin/fetch_backbone_records.py ${S}_taxonomy_nr.txt sources/genbank_backbone/records.fasta.gz
# GenBank records the backbone lacks, merged with it (the last release is the reference)
U=sources/genbank_update_20261113
python3 bin/genbank_update.py ${S}_taxonomy_nr.txt releases/nemits2_v0.1.3 2026/03/17 $U curation/genbank_named.tsv
python3 bin/genbank_update.py merge $S sources/genbank_backbone/records.fasta.gz $U work/its2
python3 bin/its2.py work/its2/seqs.fasta work/its2/records.fasta.gz work/its2/its2_cut
python3 bin/informal_accessions.py work/its2/taxonomy.txt > work/its2/informal.txt   # never a partial's relative
python3 bin/its2_partials.py work/its2/its2_cut work/its2/its2_refs 8 work/its2/informal.txt
python3 bin/ncbi_names.py ${S}_taxonomy_nr.txt curation/ncbi_names.tsv   # current NCBI names; refetch before a release
python3 bin/build.py work/its2/its2 work/its2/taxonomy.txt nemabiome_20260317 work/build --workers 8 \
    --layer curation/long_read_its2.tsv --lineages curation/lineages.tsv
python3 bin/species_coherence.py work/build > work/build/species_coherence.tsv   # screen before a release
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
