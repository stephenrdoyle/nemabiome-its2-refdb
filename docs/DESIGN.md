# nemabiome-its2-refdb: design

The ITS2 counterpart of vespa-18s-refdb: a versioned reference database
for amphid, built from the nemabiome ITS2 database rather than a primer
scan of a general one. Decisions so far, with the measurements behind
them.

## 1. Fixed inputs

- **Marker and primers.** NC1 `ACGTCTGGTTCAGGGTTGTT` (3' end of 5.8S)
  and NC2 `TTAGTTTCTTTTCCTCCGCT` (28S), from Gasser et al. 1993 (Nucleic
  Acids Res 21:2525). The sequences match nemabiome's own protocol
  page, which lists each with its Illumina adapter and 0-3 N spacers.
  Avramenko et al. 2015 report a 311-331 bp product for strongylids.
- **Reads.** Illumina 2x150 (also 2x250, nemabiome's usual MiSeq run)
  and nanopore, both on NC1/NC2. Nanopore here is not the NC5/NC2
  ITS1-5.8S-ITS2 long amplicon.
- **Scope.** All nematodes in the backbone (one phylum, 2 classes, 9
  orders, 101 families, 426 genera).

## 2. Backbone

The nemabiome ITS2 database (Workentine et al. 2020, BMC Genetics
21:74), built by markerDB: GenBank nematode records, ITS2 cut out with
an HMM, taxonomy from NCBI. Zenodo holds only v1.0.0 (2019, 8,630
non-redundant). The live web app has version 20260317 (17 March 2026,
15,108 non-redundant). It holds all but 62 v1.0.0 accessions, with
identical ITS2 for all but 2, and adds 6,540, taking genera from 325 to
426 and species from 1,429 to 1,897. 20260317 is the backbone;
`sources/README.md` says how to download it. markerDB drops "sp."
names, so every species label is a name. The only labels that do not
name one species are hybrids ("A x B"), cf./aff., sensu lato and
"environmental sample", 12 in all; they do not vote at species.

## 3. The reference unit

### 3.1 The HMM's boundaries do not match the reads

A read pair from NC1/NC2 covers 5.8S tail + ITS2 + 28S stretch, primer
to primer. For a sample of source records that run primer to primer,
the backbone's ITS2:

- starts 3 bases before NC1 ends in most records (inside the primer),
  but about 20 bases earlier in ~270, and 79 earlier in one;
- ends ~47 bases before NC2 where the record runs to NC2, but where the
  record stops sooner it often includes the 28S stretch and part of NC2.
  In ~1,700 records it does.

Whether a read and a reference are compared over the same bases would
change from record to record. All three options considered:

| | |
|---|---|
| (a) keep only records spanning NC1..NC2 | exact match to the reads, but only ~10% of records in a sample of 200 (both sites in 15) |
| (b) add a borrowed 28S stretch to ITS2-only records | keeps everything, but the 28S stretch has 20 variants among 109 primer-to-primer records, 1-6 bases apart: up to 2% invented differences, against a 0.3% species limit |
| (c) **ITS2 only, cut at fixed anchors; amphid scores reads over the reference's span** | chosen |

### 3.2 The anchors (`bin/its2.py`)

The full GenBank record behind every backbone accession is fetched
(`bin/fetch_backbone_records.py`) and the ITS2 re-cut:

- **start**: the 3' end of NC1, within 3 edits. Otherwise within 6
  edits but only within 30 bases of the backbone's start. 3,105 v1.0.0
  records carry 5.8S but differ from NC1 by 4-6 bases, nearly all
  outside the Strongylida.
- **end**: the start of the 28S stretch, matched to its commonest
  variant (`S28`, 44 bases) within 8 edits, else 12 near the backbone's
  end.
- **no anchor**: where the record starts or ends exactly at the
  backbone's boundary, the submitter trimmed there and the end is
  kept as `record`. Where flanking sequence exists but no anchor fits,
  the backbone's boundary is kept (`hmm`).
- **guard**: a cut more than 30 bases longer than the backbone's ITS2
  falls back to the backbone's boundaries (LT159837 came out 1,187
  bases long, the anchor found in the wrong place).

Of 15,108 records: 4,832 anchored at both ends, 6,268 at the start
only, 1,539 at the end only, 2,469 at neither. Both-anchored ITS2s are
94-692 bases (median 230). 1,353 records run on to NC2; their 28S
stretch (285 variants, 40-49 bases) is kept in `its2_cut.tsv` and never
added to other records.

Strongylids fit the primers: in Trichostrongylidae, Strongylidae and
Ancylostomatidae 85-90% of records holding the NC1 site match it
exactly. Most of the rest begin partway through the primer, not true
mismatches.

## 4. Partial records (`bin/its2_partials.py`)

Each record not anchored at both ends is mapped against those that are
(minimap2; see 4.1). The best relative, at >= 80% identity, stands in
for the missing anchor:

- **trim**: bases past the relative's ITS2 are cut off, if the
  alignment reaches the relative's end (within 3 bases) and at least 8
  hang over. 2,406 records lose a 28S stretch stub at the end (median
  23 bases). 443 lose 5.8S or primer remnants at the start (median 48).
  The 5.8S bases would cost a read a mismatch, as reads do not have
  them. Checked by hand: the kept sequence starts where the relative's
  does.
- **grade** by coverage of the relative's ITS2: `near_complete` >= 98%,
  `partial_90`, `partial_75`, `partial` < 75%; `unmeasured` with no
  relative at 80%.

80% is low on purpose. Sister species within a genus (*Labiosimplex*)
are often 85-90% identical over ITS2, yet still align end to end, and
that is all a boundary estimate needs.

| Grade | Records | Species whose best record it is |
|---|---:|---:|
| complete | 4,832 | 871 |
| near_complete | 4,889 | 362 |
| partial_90 | 1,333 | 76 |
| partial_75 | 1,009 | 38 |
| partial | 1,027 | 33 |
| unmeasured | 2,018 | 517 |

Unmeasured records are mostly plant-parasitic and free-living families
with few complete ITS2s: Longidoridae, Steinernematidae,
Aphelenchoididae, Hoplolaimidae, Rhabdiasidae. Their median length,
296, suggests many are whole ITS2s. Partials are kept, graded, not
dropped: a species missing from the database still sends reads, which
land on the nearest remaining species.

### 4.1 minimap2's repeat filter

By default minimap2 ignores minimizers shared by too many targets.
Thousands of near-identical strongylid ITS2s push the common ones over
that limit. Partials then found no relative at all, *Strongylus
vulgaris* among them. Every minimap2 call here uses `-f 1000000`. (The
same effect dropped pairs in vespa-18s-refdb's neighbour check.)

## 5. Collapse and release (`bin/build.py`)

- **Merge**: identical ITS2s become one record. Each rank takes the
  single informative name its members agree on; two names stop it
  there. *B. malayi* and *B. pahangi* share an ITS2, so their record
  stops at *Brugia*.
- **Conflicts**: disagreement at order or above is treated as a
  mislabel: a 2/3 majority of at least two wins, otherwise the record
  merges up and is listed. There are none in 20260317: identical ITS2s
  never cross an order.
- **Grade** of a record = its best member's.
- **Windows**: a 2xR read pair with primers trimmed sees the first R-20
  bases of the ITS2 and the last R-67 (R-20 less the 28S stretch):
  130 + 83 for 2x150, 230 + 183 for 2x250. A record no longer than the
  two windows is seen whole.
- **Margins** (`species_margins.tsv`): for each species-level record,
  the edits to the nearest record of another species. Each record is
  placed inside the other, with its overhangs free, and the better
  direction taken, so a partial is compared only over the bases it
  has, as amphid will score it. The first build measured margins over
  both whole sequences instead. That counted a partial's missing bases
  as differences, and put 99% of records at species.
- **Ties**: a species-level record at margin 0 sits exactly inside
  another species' record over its span. A read scores the same
  against both, so the resolution tables stop it at the rank the two
  share. 732 records read whole, 745 through 2x150 windows.
- **Windows never resolve deeper than whole**: the build stops if any
  would. The first run did stop. A record already at genus was named to
  species through a window group, because a blank species was read as
  "no vote" rather than "these members disagree". Fixed in
  `merge_records`.

### 5.1 The first build (scratch, 2026-11-09, curated)

`results/dev_build_summary.txt`. 12,980 records, 1,838 species with a
record of their own (1,831 through 2x150 windows), 424 genera. Read
whole, 94% of records resolve to species and 6% to genus; through 2x150
windows, 92% and 7%.

| Family | Records | Species | To species, whole | 2x150 |
|---|---:|---:|---:|---:|
| Strongylidae (incl. cyathostomins) | 1,343 | 265 | 100% | 99% |
| Trichostrongylidae | 1,153 | 66 | 93% | 93% |
| Ancylostomatidae | 288 | 17 | 97% | 97% |
| Metastrongylidae | 454 | 47 | 100% | 99% |
| Onchocercidae | 385 | 24 | 97% | 93% |
| Trichuridae | 339 | 19 | 99% | 94% |
| Ascarididae | 121 | 30 | 77% | 67% |
| Heteroderidae | 1,148 | 109 | 87% | 87% |

"Resolves to species" means the record is distinct from every other
species' record. It does not mean reads are named correctly. Margins
are thin: of 12,855 species-level records, 705 are tied with another
species, 1,019 are one edit away, 836 two and 1,628 three to five. Over
a 240-base ITS2, one edit is 0.4%. Intraspecific ITS2 variation is
often that large, and in vespa-18s-refdb a leave-one-out test named
27-46% of reads as a sister species once identity limits were lifted.
The same test is due here (§7).

### 5.2 Cross-family ties: mislabels, curated

The uncurated build had 51 species-level records sitting exactly inside
a record of another family (`results/cross_family_ties.tsv`). Every
record on either side, 59 in all, was BLASTed against nt (megablast).
Each record's label was then checked against the other GenBank entries
at >= 98% identity over >= 90% of it
(`results/cross_family_ties_blast.tsv`). In every case one side's label
was contradicted by all the other near-identical entries. The 14
mislabelled source records are dropped in `curation/overrides.tsv`. A
15th was dropped on a second pass, after the first 14 drops exposed it:
GQ258655, "*Steinernema siamkayai*", whose ITS2 is identical to three
*Heterorhabditis* records.

- *Contracaecum muraenesoxi* MW420927, MW412588, EU826125: ITS2 of
  *Hysterothylacium amoyense* / *H. zhoushanense*. 15 partial
  *H. amoyense* records sat inside the first.
- *Heterorhabditis* / *Steinernema* / *Oscheius* / *Acrobeloides*
  (entomopathogenic and free-living nematodes): 8 records, including an
  apparent label swap within one submission (MN704735 / MN704744).
- *Oscheius tipulae* PX455932 = *Steinernema glaseri*; *Pratylenchus
  crenatus* PQ403988 = *Punctodera*; "*Caenorhabditis elegans*"
  EU131007 (NCBI: UNVERIFIED_ORG) = *Pratylenchus speijeri*.
- "*Cercopithifilaria*" JF501396 = *Crenosoma vulpis*, a dog lungworm.
  Probably the same dog's sample.
- *Parabronema skrjabini* MK760918 inside a *Marshallagia marshalli*
  record.
- ***Oesophagostomum radiatum* MG685662 = *Mecistocirrus digitatus*.**
  *O. radiatum* is on vespa-18s-refdb's watch list; its 5 other ITS2
  records are unaffected.

None remain after curation. One record still stops at order and is not
a tie: two *Dujardinascaris gigantea* and one *Ortleppascaris sinensis*,
both crocodilian ascaridoids that NCBI puts in different families. It
is either a misidentification or shared ITS2; it is left, merged.

### 5.3 Neighbour check

Exact ties miss a mislabel that differs from its true species by a base
or two. As in vespa-18s-refdb's neighbour check, every species-level
record whose nearest other-species record (over the shorter span) was
>= 97% identical across families, or >= 99% across genera, was BLASTed
with its neighbour: 209 pairs, 219 records
(`results/neighbour_check_blast.tsv`). Seven more source records are
dropped, each contradicted by 10-49 other GenBank entries:

- *Cactodera estonica* AY692355 = *Heterodera*. This one record made 66
  of the 209 pairs.
- *Anisakis typica* PX761267 = *Hysterothylacium*.
- *Ascaridia galli* MH392317 = *Heterorhabditis*.
- *Contracaecum muraenesoxi* EU828749, MW418099, FJ234209 =
  *Hysterothylacium amoyense*. With the three in §5.2, every *C.
  muraenesoxi* ITS2 in the backbone is dropped and the species leaves
  the database.
- *Acrobeloides saeedi* MK935151 = *Heterorhabditis indica*.

97 pairs remain (78 across genera >= 99%, 19 across families >= 97%):

- **The same species under two genus names**: *Ostertagia trifurcata* /
  *Teladorsagia circumcincta* (male morph); *Phocanema* /
  *Pseudoterranova*; *Labiosimplex* / *Labiostrongylus longispicularis*;
  *Parapoteriostomum* / *Cylicodontophorus mettami*; *Gracilacus* /
  *Paratylenchus paralatescens*; *Betulodera* / *Cactodera* /
  *Heterodera betulae*; *Japanema* / *Sheraphelenchus entomophagus*.
  To be merged by relabels once the names are chosen.
- **Unclear**:
  - *Rostellascaris spinicaudatum* vs *Toxocara cati* (2-3%);
  - *Oesophagostomum aculeatum* LC428824 vs *Ternidens deminutus*
    (watch list; 99%);
  - *Rugopharynx alpha* vs *Papillostrongylus barbatus*;
  - *Steinernema khoisanae* vs *Oscheius chongmingensis*;
  - *Ektaphelenchoides* vs *Bursaphelenchus*.

  Left in place.
- **Cyathostomins**: *Cyathostomum catinatum* and *Cylicostephanus goldi*
  records are 99.4-99.6% identical. HiFi reads of both species (§5.4)
  show the genuine variants of the two are 4-8 edits apart. The one
  *C. goldi* record behind the flag (PP923309) is 2 edits from the
  nearest *C. catinatum* read variant and 3 from the nearest *C. goldi*
  one: too close to call, so it is left.

### 5.4 Genomes and HiFi reads

`bin/its2_scan.py` finds the ITS2 in an assembly or a read set. It can't
work as a primer scan does for V4: NC1 is 4-6 bases off in many orders,
and no 5.8S k-mer is shared by all of them (Trichinellida lack
GAAGAACGCAGC). So the rDNA is located by mapping to a panel: up to two
5.8S-ITS2-28S stretches from each of 258 genera, taken from records
anchored at both ends. Each hit is then cut at the NC1 and 28S anchors
used for the backbone (§3.2). Assemblies are cut into 20 kb windows;
reads are taken whole.

Four faults found and fixed on the way:
- **minimap2's chaining thresholds**: at the defaults (`-m/-s 40`) it
  missed *Spirocerca*, which shares only the 5.8S tail and 28S head with
  its nearest panel relative.
- **Assembly gaps**: the matcher read N as any base, so runs of N
  "matched" both anchors (ten 50-base ITS2s of pure N in *Dracunculus*).
- **Soft-masked assemblies**: repeats in lower case counted as
  mismatches throughout.
- **Reads starting inside the ITS2**: these took a 6-edit NC1 look-alike
  and gave truncated ITS2s (11% of *Bunostomum* reads). A loose NC1 is
  now accepted only when the alignment starts at or before NC1 on the
  panel record, and nothing shorter than 80 bases is reported.

**Gap species** (watch-list nematodes with no ITS2 in the backbone). The
assemblies that gave V4 were scanned:

| Species | Assembly | ITS2 copies (distinct) | Edits between variants | Nearest backbone record |
|---|---|---:|---|---|
| *Loa loa* | GCA_053739935 | 3 (2) | 1 | *Dirofilaria* 64% |
| *Acanthocheilonema viteae* | GCA_046563165 | 13 (9) | 6 median, 15 max | *A. reconditum* 86-88% |
| *Litomosoides sigmodontis* | GCA_963070105 | 5 (5) | 27 median, 37 max | none above 56% |
| *Spirocerca lupi* | GCA_965124195 (Tree of Life) | 24 (9) | 17 median, 46 max | filarioids ~61% |
| *Trichuris incognita* | STPHxBioinformatics | 4 (4) | 28 median, 31 max | *T. trichiura* 96-100% |
| *Dracunculus medinensis* | WBPS19 | 0 | - | - |

- *Dracunculus*: NC1 sits about 6 bases off and NC2 about 6, but no 28S
  stretch fits. Its rDNA is too far from any panel genus; there is no
  Dracunculoidea or Camallanina in the backbone.
- One *T. incognita* copy is identical to a GenBank "*T. trichiura*"
  ITS2. The human *Trichuris* records may include *T. incognita*.
- None of these is added yet. Copies from older short-read assemblies
  (*A. viteae*, *Litomosoides*, *T. incognita*) may carry assembly
  errors; the *Spirocerca* Tree of Life HiFi assembly is the most
  trustworthy.

**HiFi read sets** (the ten runs used for V4 in vespa-18s-refdb, first
200,000 reads each). A variant is kept if it is >= 2% of the run's ITS2
reads and seen in >= 5 reads:

| Species | ITS2 reads | Variants | Edits between them | Main variants vs backbone |
|---|---:|---:|---|---|
| *Bunostomum trigonocephalum* | 2,449 | 3 | 1-2 | same species, 0-2 edits |
| *Cyathostomum catinatum* | 1,402 | 7 | 1-4 | same species, 0 (one 2) |
| *Cylicostephanus goldi* | 3,789 | 3 | 2-8 | same species, 0-1 |
| *Cylicostephanus longibursatus* | 1,442 | 14 | 1-12 | same species, 0-1 (one 7 from anything) |
| *Cylicocyclus insigne* | 496 | 6 | 1-6 within each of two groups ~20 apart | *C. insigne* and *C. ashworthi*, 0 |
| *Cylicocyclus nassatus* | 933 | 2 | 2 | same species, 0 |
| *Coronocyclus labiatus* | 183 | 1 | - | same species, 0 |
| *Ancylostoma duodenale* (subreads) | 95 | 5 | 2-8 | same species, 0-4 |
| *Angiostrongylus malaysiensis* | 4,105 | 8 | 2-8, plus a group 22-30 away | same species, 0-4 |
| *Crenosoma striatum* | 781 | 8 | 1-4, plus 2 minor variants 11-24 away | same species, 0-2 |

What the reads show:

- **Labels hold.** Every run's main variants match backbone records of
  the same species, usually exactly. That is independent confirmation
  of the GenBank labels for nine of the ten species; the tenth is
  *C. insigne*, below.
- **Variation within one worm is as large as between species.** Variants
  1-12 edits apart are routine; *C. goldi* and *C. catinatum* variants
  are 4-8 apart, *C. longibursatus* 7-9 from both. Among cyathostomins
  no identity or edit limit separates species. A read can be named to
  species only when the database holds its exact variant. Many variants
  already are in GenBank: most read variants match a record exactly.
- ***Cylicocyclus insigne* / *C. ashworthi*.** ERR11658755 is one adult
  worm (Sallé et al., PRJEB63274; INRAE, identified by ITS2 and COX1).
  It carries two ITS2 types ~20 edits apart: a 324-base type matching
  *C. insigne* records (34% of reads) and a 318-base type matching
  *C. ashworthi* records exactly (43%). Either *C. insigne* carries an
  "*ashworthi*" copy type, or the worm is a hybrid or misidentified;
  its identification by ITS2 makes this circular. Reads matching the
  *ashworthi* type may come from *C. insigne*. Left as is, flagged.
- **Divergent minor copy types.** *A. malaysiensis* and *C. striatum*
  carry a small group of variants 11-30 edits from the main ones. They
  are either a second rDNA copy type or a second worm in the sample.

Not yet added to the database. A read-variant layer, as in
vespa-18s-refdb, would put the exact within-worm variants in, which is
what species calls depend on (§7).

## 6. Release format

As vespa-18s-refdb (`refdb.fasta`, `taxonomy.tsv`, `members.tsv`,
`conflicts.tsv`, `primer_fit.tsv`, resolution tables, `anchor_lookalikes.tsv`,
`manifest.json`, `summary.txt`), with these differences:

- taxonomy columns are NCBI's eight ranks (superkingdom ... species),
  not PR2's nine positional levels;
- `grade` and `coverage` in `taxonomy.tsv` and the resolution tables;
  `members.tsv` records the bases trimmed at each end;
- `species_margins.tsv` is written by the build itself, with the nearest
  record per geometry;
- `primer_fit.tsv`: `exact`, `near` (<= 2 edits), `poor`, `one_site`
  (the source record holds only one primer site), `unknown`. 87% of
  records are `one_site` or `unknown`: GenBank ITS2 records seldom run
  primer to primer;
- `manifest.json` adds `score_over`, the rule the windows and margins
  assume.

amphid needs changes before it can use this database:
[feedback/its2_notes_for_amphid.md](feedback/its2_notes_for_amphid.md).

## 7. Open questions

1. **Ties within a family** (27 records stop at family, read whole):
   check as for §5.2? Cross-genus ties are more often real than
   cross-family ones.
2. **Leave-one-out simulation** for ITS2, as in vespa-18s-refdb: how
   often a read from a species missing from the database is named as a
   sister species, at amphid's limits and at ITS2-appropriate ones.
3. **Intraspecific variation**: answered in part by §5.4. Within one
   worm, variants differ by 1-12 edits, as much as some species do. A
   species call needs the read's exact variant in the database, not an
   identity limit. Next: a read-variant layer from the HiFi runs and the
   *Spirocerca* assembly, and the same scan of *Haemonchus*,
   *Teladorsagia* and *Trichostrongylus* genomes.
4. **Contained records of the same species.** Partials identical to
   part of a same-species complete record are separate records. That is
   harmless for calls, but it inflates record counts. Fold them into
   the complete record?
5. **Unmeasured records** (2,018): a lower-identity relative, or an HMM
   for the 28S stretch, could anchor more of them.
6. **ITS1-5.8S-ITS2** for NC5/NC2 nanopore runs, should those be used
   later: nemabiome has a separate database (v0.0.1).
