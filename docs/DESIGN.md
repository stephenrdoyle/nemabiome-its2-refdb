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
- **part**: a record that starts with the tail of NC1 (the last >= 5
  bases exact, or >= 8 with one mismatch) or ends with the head of the
  28S stretch has that remnant cut off (`part`). Before this, 155
  records started inside NC1 and 494 ended inside 28S, and a read was
  charged for bases it does not carry. A remnant with an indel is not
  caught: AM410070 (*H. placei*) keeps GGGTGTT at its start.
- **guard**: a cut more than 30 bases longer than the backbone's ITS2
  falls back to the backbone's boundaries (LT159837 came out 1,187
  bases long, the anchor found in the wrong place).
- **strand**: some records are deposited reverse-complemented, and the
  anchors cannot tell, as a reversed record may still fit a loose anchor.
  Records with a tight NC1 or 28S anchor (`primer`) set the strand;
  every other cut is mapped against them (minimap2) and, if its best
  match is on the minus strand, turned and cut again. 36 records turned,
  among them MK087176 (*H. contortus*: 40 edits from its own species
  before, 2 after) and AM410070. LN626259 (*Dirofilaria*) maps equally
  well both ways and is left as deposited.

Of 15,108 records (counts from before the strand and `part` steps): 4,832 anchored at both ends, 6,268 at the start
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
  share. 343 records read whole, 473 through 2x150 windows (the current
  v0.1.1; 410 and 522 in v0.1.0, 732 and 745 before curation
  and folding).
- **`species_status`** (`taxonomy.tsv`): `named` if the record reaches
  species; `shared` if it stops above species because its members name
  two or more species (95 records: real shared ITS2s, such as *Brugia*);
  `unnamed` if no member names a species at all (32: hybrids, cf.,
  environmental samples). The classifier's gap counts a shared record as
  another species - a read matching it could be either - but ignores an
  unnamed one, which says nothing about which species a read is from.
  Before this, three hybrid records cut species calls in a cattle set
  from 87% to 35%.
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
The same test is in §6.

v0.1.0, after the curation of §5.2-5.3, folding (§8, item 4) and the
long-read layer (§5.5), has 11,956 records, 1,840 species with a record
of their own and 426 genera; 96% of records resolve to species read
whole, 94% through 2x150 windows. v0.1.1 adds
the strand and `part` steps (§3.2) and folds weak records (§8, item 9):
11,146 records, the same species and genera, 96% to species read whole
and 94% through 2x150 windows (`results/dev_build_summary.txt`).

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

- **The same species under two names**: *Ostertagia trifurcata* /
  *Teladorsagia circumcincta* and *O. lyrata* / *O. ostertagi* (minor
  male morphs; relabelled in `curation/overrides.tsv`); *Parapoteriostomum* /
  *Cylicodontophorus mettami* and *Betulodera* / *Cactodera betulae*
  (settled by current NCBI names, §5.6). NCBI keeps four look-alike pairs
  as separate species, and so does this database: *Phocanema* /
  *Pseudoterranova*, *Labiosimplex* / *Labiostrongylus longispicularis*,
  *Gracilacus* / *Paratylenchus paralatescens*, *Japanema* /
  *Sheraphelenchus*.
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

**The *Haemonchus contortus* / *H. placei* shared variant.** One ITS2 was
deposited twice as *H. contortus* (KX829170, MH481602) and twice as
*H. placei*. It is the most common ASV in a cattle larval set
(PRJNA1524548, 23% of reads), which therefore stopped at genus. Two
*H. contortus* assemblies hold 13 ITS2 copies between them: MHco3(ISE)
(WBPS19) and GCA_041937105.1. Every copy is 5-9 edits from that variant
and matches *H. contortus* records exactly; every copy is 6-10 edits
from *H. placei* records. *H. contortus* rDNA does not carry the variant.
Both "*H. contortus*" records came from cattle, *H. placei*'s host, so
they are probably *H. placei* or hybrids. Both are dropped, and the
variant is called *H. placei*. Four more "*H. contortus*" records from
the same cattle survey as MH481602 (MH481579, MH481580, MH481588,
MH481596) are 1 edit from *H. placei* records and 5-6 from every
*H. contortus* genome copy; they are dropped too. The *H. placei* assembly (WBPS19
PRJEB509, short reads) gave no ITS2, nor did McMaster or GCA_007637855.2.

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
| *Bunostomum trigonocephalum* | 2,151 | 3 | 1-2 | same species, 0-2 edits |
| *Cyathostomum catinatum* | 1,411 | 7 | 1-4 | same species, 0 (one 2) |
| *Cylicostephanus goldi* | 3,789\* | 3 | 2-8 | same species, 0-1 |
| *Cylicostephanus longibursatus* | 1,438 | 14 | 1-12 | same species, 0-1 (one 7 from anything) |
| *Cylicocyclus insigne* | 1,620 | 7 | 1-7 within each of two groups 18-23 apart | *C. insigne* and *C. ashworthi*, 0 |
| *Cylicocyclus nassatus* | 719 | 2 | 2 | same species, 0 |
| *Coronocyclus labiatus* | 182 | 1 | - | same species, 0 |
| *Ancylostoma duodenale* (subreads) | 798 | 5 | 2-8 | same species, 0-4 |
| *Angiostrongylus malaysiensis* | 4,103 | 8 | 2-8, plus a group 22-30 away | same species, 0-4 |
| *Crenosoma striatum* | 730 | 8 | 1-4, plus 2 minor variants 11-24 away | same species, 0-2 |

\* From the first scan; the rescan with the fixed scanner got only 67,323
reads before EBI's download stopped (1,300 ITS2 reads, the same three
variants). *C. nassatus* and *C. striatum* stopped at 155,204 and 189,087
reads.

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
  It carries two ITS2 types 18-23 edits apart: a 324-base type matching
  *C. insigne* records and a 318-base type matching
  *C. ashworthi* records exactly. Either *C. insigne* carries an
  "*ashworthi*" copy type, or the worm is a hybrid or misidentified;
  its identification by ITS2 makes this circular. Reads matching the
  *ashworthi* type may come from *C. insigne*. Left as is, flagged.
- **Divergent minor copy types.** *A. malaysiensis* and *C. striatum*
  carry a small group of variants 11-30 edits from the main ones. They
  are either a second rDNA copy type or a second worm in the sample.

These copies, and the assemblies' below, are now in the database as a
layer (§5.5).

### 5.5 The long-read layer

`curation/long_read_its2.tsv`, made by `bin/genome_layer.py` from the
scans listed in `curation/long_read_sources.tsv` and merged by
`bin/build.py --layer` like any other source. Every copy is cut at the
backbone's anchors, so it collapses with identical backbone records.

- **Assemblies** (31: 24 long-read assemblies of 16 species, the two
  *H. contortus* assemblies of §5.3, and the five gap-species assemblies
  of §5.4 that gave an ITS2): every distinct copy is kept, except that
  in an array of >= 20 copies a variant must be seen twice (one copy
  in a big array may be an assembly error).
- **HiFi reads** (the ten runs above): a variant is kept if it is >= 2%
  of the run's ITS2 reads and seen in >= 5.
- **Labels** come from the backbone's records of the species (else its
  genus), else NCBI Taxonomy.

174 copies: 51 match a backbone record, 97 are new records of 25
species, and
six species enter the database: *Loa loa*, *Acanthocheilonema viteae*,
*Litomosoides sigmodontis*, *Spirocerca lupi*, *Trichuris incognita*,
*Heligmosomoides bakeri*. The two pre-release Tree of Life
*Angiostrongylus* assemblies are not used.

Of the 51 matches, 47 are to records of the same species - more
confirmation of the labels. The other four stop at genus: three
*C. insigne* read variants identical to *C. ashworthi* records (the
shared type of §5.4), and a *Parascaris univalens* copy identical to a
record shared with *P. equorum*.

### 5.6 Names: current NCBI Taxonomy

markerDB named each record when the backbone was built; NCBI has renamed
taxa since. Rather than choose names pair by pair, the rule is: **use
NCBI Taxonomy's current name**. `bin/ncbi_names.py` fetches the current
order, family, genus and species of every backbone taxid (1,915) into
`curation/ncbi_names.tsv`, and `bin/build.py --names` applies it before
`curation/overrides.tsv` (which keeps the last word). A rank is renamed
only where the record names it, and names that are not one species
(hybrids, cf., environmental) are left alone. Each renamed source record
says so in `members.tsv`.

On 2026-11-10, 64 source records (20 species names, 14 genera, 2
families) took a new name, among them:

- three synonym pairs that had split one species in two: *Parapoteriostomum
  mettami* -> *Cylicodontophorus mettami*, *Betulodera betulae* ->
  *Cactodera betulae*, *Tylenchorhynchus ventrosignatus* -> *Bitylenchus
  ventrosignatus* (named species 1,872 -> 1,869). The *mettami* split had
  stopped the k-mer classifiers at family on 5% of equine mock reads (§6.3);
- genus moves such as *Phasmarhabditis* -> *Pellioditis* (4 species),
  *Geocenamus*/*Merlinius* -> *Scutylenchus*, *Heterakis spumosa* ->
  *Ganguleterakis spumosa*, *Aphelenchoides avenae* -> *Aphelenchus
  avenae* (Aphelenchidae).

Resolution is unchanged except for one record: NCBI moved *Oscheius
microvilli* to *Caenorhabditis*, but its ITS2 is identical to an
*Oscheius myriophilus* record, so the shared record now stops at
Rhabditidae. ITS2 sides with *Oscheius*; the rule follows NCBI, and the
case is left as a flag. The evaluation in §6 was run on the build before
this step; renaming changes labels, not sequences or ties.

## 6. Naming reads: the classifier and its tests

amphid cannot use this database yet (§7), so a small stand-in scores
reads the way the amphid notes ask: `bin/its2_classify.py` for raw read
pairs, `bin/classify_asvs.py` for DADA2 ASVs. Both share the scorer.

- **Span scoring.** Each reference is placed whole inside the read, the
  read's extra bases free. 2x150 pairs that do not overlap are scored
  as two windows, against the start and the end of each reference.
  Every reference base from its own end to where the window stops is
  charged; read bases before the reference's start are free only up to
  30 (primer offset and flank). The second window is cut at the 28S
  anchor first (`trim_28s`, else 47 bases): the stretch runs 40-49
  bases, so a fixed cut left a third of records' last bases nothing to
  meet. v0.1.0's scorer gave up to 10 of a record's first bases free
  (`SLACK`), which let a partial *H. bedfordi* record (§6.2) tie reads
  of *H. placei*.
- **Strand.** An ASV whose best minimap2 hit is on the minus strand is
  turned round before scoring (none of the 810 real-data ASVs; one
  mock ASV, §6.3).
- **Every candidate scored.** minimap2 proposes; every reference it
  proposes is scored exactly. v0.1.0's figures came from the top 40 by
  matching bases, which hid the nearest other species from reads of
  species with hundreds of records: that neighbour is often a partial
  with fewer matches (issue #1). The figures below are uncapped, on the
  v0.1.1 (§5.1).
- **Ranking by identity**, not raw edits, over at least 100 bases:
  ranked by edits, a 77-base partial 32 edits off beat full-length
  references.
- **Ties**: every reference at the best identity whose span is within
  5 bases of the longest such span. A partial matching only over a
  shorter span does not tie with a full-length match.
- **The gap rule**: species when >= 99% identical and the nearest
  reference of any other species is >= 4 edits further off, or when the
  match is exact and it is >= 2 edits off (sister species share
  variants: an exact match with another species 1 edit off is right
  60% of the time, §6.1); genus at >= 97%, family at >= 95%. `shared` records count as another species,
  `unnamed` ones do not (§5).
- **What it reports**: identity, ties and gap, and two columns that say
  how close a call that stopped short came. `nearest` is the tied
  species and their identity; `next` is the nearest other species and
  how many edits further it is. A shared record is shown with the
  species its members name (*Cooperia* shared: punctata / spatulata).

### 6.1 Leave-one-out

As in vespa-18s-refdb. `bin/its2_leave_one_out.py`: up to 5
full-length species-level records per species, 3,195 queries from 1,224
species, read whole and through 2x150 windows, as perfect reads.

- **Missing**: every record holding a member of the query's species is
  removed. Any species named is wrong; genus or above is right.
- **Present**: only the query's own record is removed, for species with
  another record.

Reported as recall (right calls / queries), precision (right calls /
calls made at that rank or deeper; a species call is right at genus if
its genus is) and over-classification (a species named when the species
is missing; at genus, a genus named when the genus is missing too).
Specificity is not used: with 1,224 species every query is a negative
for all but one, and every method scores above 99.9%.
(`bin/loo_metrics.py`, `results/leave_one_out/metrics.txt`.)

| Rule (whole) | Species recall | Species precision | Genus recall | Genus precision | Over-classification, species | Over-classification, genus |
|---|---:|---:|---:|---:|---:|---:|
| amphid's V4 limits (99.7% / 99.5% / 95%) | 17.9% | 93.7% | 49.7% | 99.3% | 2.0% | 1.1% |
| exact: species only at 0 edits | 10.7% | 92.2% | 52.6% | 99.3% | 1.6% | 1.1% |
| gap (>= 99% and 4 edits clear, or exact and 2 clear) | 53.4% | 98.5% | 89.2% | 99.3% | 6.4% | 5.4% |

At 2x150 the gap rule gives 56.3%, 98.7%, 90.2%, 99.4%, 7.0% and 5.0%
(`results/leave_one_out/`). Before exact matches were named at 2 edits
clear, it gave 52.3%, 98.6% and 6.3% whole, 54.0%, 99.0% and 6.9% at
2x150 (`bin/loo_gap_split.py`, `results/leave_one_out/gap_split.txt`). v0.1.0, scored uncapped, gave 52.6% recall
and 6.5% over-classification read whole, 56.1% and 9.0% at 2x150: the
stricter window scoring costs 2 points of recall at 2x150 and saves a
quarter of the over-classification. The V4 and exact rules seldom
name a species even when it is there. The gap rule names three to five
times as many, and is wrong less often on present species, because the
margin to the next species is what separates right from wrong, not the
identity. Its missing-species errors are mostly reads identical to a
record of another species - a shared ITS2 deposited under two names -
which no threshold can catch.

A third of present-species reads stop at genus. 27% of queries are
below 99% identity to their own species (intraspecific variation; §5.4)
and 15% have a gap of 0-3 edits (sister species whose ITS2 barely
differs). Could a rule per genus or family do better?
`bin/loo_group_rules.py` let each genus or family choose its own
identity floor (97-99%) and gap (1-10 edits), chosen on four fifths of
the species and scored on the rest. On held-out species, per-group rules
recovered 2.3-4.4 points of right calls for 0.4-0.6 more present-wrong
and 0.4-1.0 more missing-wrong. Loosening the global rule trades 9.1
points for 1.0 and 2.9 (`results/leave_one_out/group_rules.txt`). No
per-group rule beats the global one by more than fractions of a point,
from five splits of 1,224 species; one rule is kept.

### 6.2 Against DADA2

The standard nemabiome workflow is DADA2: cutadapt, filterAndTrim,
learnErrors, dada, mergePairs, removeBimeraDenovo, then assignTaxonomy
(minBoot 80) or, as nemabiome.ca now runs it, DECIPHER IdTaxa (threshold
60, 100 bootstraps). `bin/dada2_nemabiome.R`, `bin/assign_taxonomy.R`,
`bin/idtaxa.R` (dada2 1.38.0, DECIPHER 3.6.0). 2x150 pairs do not
overlap over a 290-330 base amplicon and are joined with
`justConcatenate`. 2x250 pairs are joined hybrid: merged where they
overlap by 12 bases, concatenated where they do not (an ITS2 over about
400 bases, below).

**Leave-one-out, same queries** (`bin/loo_dada2.py`; assignTaxonomy and
IdTaxa retrain per fold, so 10 folds rather than one species at a time;
`results/leave_one_out/`):

| Classifier, this database (whole) | Species recall | Species precision | Genus recall | Genus precision | Over-classification, species | Over-classification, genus |
|---|---:|---:|---:|---:|---:|---:|
| assignTaxonomy | 85.4% | 94.6% | 97.4% | 99.2% | 43.9% | 20.7% |
| IdTaxa | 76.9% | 96.7% | 93.5% | 99.4% | 22.7% | 8.0% |
| gap rule | 53.4% | 98.5% | 89.2% | 99.3% | 6.4% | 5.4% |

At 2x150: assignTaxonomy 86.4% recall, 95.3% precision, 39.4%
over-classification; IdTaxa 77.6%, 96.9%, 22.1%; the gap rule 56.3%,
98.7%, 7.0%. The three sit on one trade-off line. At genus they are
equally precise; at species they differ in how often they name a
species that is not there. assignTaxonomy never declines: a read from a
missing species gets a species name 44% of the time. IdTaxa halves
that. The gap rule over-classifies least, by 3-7 times, and names
fewest. If 20% of a sample's reads come from species missing from the
database, the share of species calls that are right is about 84%
(assignTaxonomy), 90% (IdTaxa) and 96% (gap rule). Of missing-species
reads whose genus is still in the database, the gap rule names the
genus for 28% (precision 93%), IdTaxa 42% (93%), assignTaxonomy 71%
(93%); most of the gap rule's rest stay unassigned. Reads left
unassigned by IdTaxa and the gap rule are not junk: they are reads whose
nearest reference is too far for a family (< 95%).

**Where the k-mer classifiers' extra species calls come from**
(`bin/loo_kmer_vs_gap.py`, `results/leave_one_out/kmer_vs_gap.txt`;
whole reads, points of all present queries, right / wrong species calls
on the reads the gap rule stops short on):

| Why the gap rule stops short | Share | assignTaxonomy | IdTaxa |
|---|---:|---:|---:|
| own species best, but < 99% | 22.9 | 21.1 / 0.2 | 18.1 / 0.1 |
| own species >= 99%, a sister 0-3 edits off | 11.6 | 8.7 / 1.0 | 6.3 / 0.2 |
| another species closest | 9.2 | 2.8 / 2.7 | 0.5 / 1.6 |
| tied with another species | 2.4 | 1.1 / 0.4 | 0.2 / 0.0 |

- **The 99% floor is most of it.** Two thirds of assignTaxonomy's
  33-point recall lead are reads from a divergent variant of a present
  species, which it names almost always correctly. Its confidence is
  relative (does this species keep winning the bootstrap?), the gap
  rule's absolute (is the read within 1% of a record?). The same
  relative judgement names missing species: of hidden queries whose
  nearest record is another species, assignTaxonomy names a species for
  36 points of all hidden queries, IdTaxa for 17.
- **One or two diagnostic sites suffice for them.** The gap rule wants 4
  edits clear of the next species (2 for an exact match); the k-mer classifiers resolve sisters
  1-3 edits apart, mostly correctly when the species is present. The
  real-data genus stops are this kind (*Cooperia* +1).
  Tuned on held-out species, the gap stays at 4 edits or more for whole
  reads and at 3 or 4 at 2x150 (`results/leave_one_out/scaled_floor.txt`,
  plain rule, lam 4).
- **Where another species is closest, they gain nothing**: assignTaxonomy
  is right as often as wrong, IdTaxa wrong more often than right.

They do not see what identity misses. Within one identity band of the
best hit, the k-mer classifiers name a species for a present species'
variant more often than for a missing species, but not by much: at
98-99%, assignTaxonomy 95% against 68% (338 right calls to 171 wrong),
IdTaxa 88% against 57%; below 98%, wrong calls outnumber right ones.
About 2 right calls per wrong one is also what a distance-scaled floor
gives the gap rule: below 99%, a species at >= 98% when the next species
is >= max(6, 3 x edits to the best) further off. Chosen on four fifths
of the species and scored on the rest, as for the per-group rules
(§6.1), it adds 7 points of recall at unchanged present-species
precision for 2.5 more points of over-classification, and is chosen
only when a wrong call is valued at 2 right ones or fewer; at 4 or 8
every fold keeps the current rule (`bin/loo_scaled_floor.py`,
`results/leave_one_out/scaled_floor.txt`). So the difference is where
each sets its bar, not information: the k-mer classifiers are tuned for
the species being in the database, the gap rule for it perhaps not
being. The current rule is kept, with no relaxed mode.

**Real data**: two published NC1/NC2 sets, sequenced 2x251.
Francis & Šlapeta (2023, Parasitology; PRJNA936713): larval cultures
from 10 cattle, 10 feral goats and 10 sheep on one New South Wales
property. Gilleard (PRJNA1524548): 82 libraries, cattle larval cultures
and faecal samples from kangaroos, pigs and rabbits, two libraries per
sample. Each run from raw reads by `bin/real_data.sh`: primer trimming
(`bin/trim_primers.py`), DADA2 denoised once and joined three ways
(merged, concatenated, hybrid), and the same reads cut to 2x150. Share
of reads named to species (`results/dada2_comparison/`, with an HTML
report per set):

| Set | Reads | assignTaxonomy x nemabiome | IdTaxa x nemabiome | assignTaxonomy x this db | IdTaxa x this db | gap rule x this db |
|---|---:|---:|---:|---:|---:|---:|
| Francis & Šlapeta 2x250, hybrid | 488,062 | 72.0% | 81.7% | 99.5% | 99.4% | 92.0% |
| merged only | 465,576 | 70.2% | 81.2% | 99.3% | 99.7% | 92.9% |
| 2x150 | 484,432 | 76.0% | 88.2% | 98.9% | 97.8% | 92.4% |
| Gilleard 2x250, hybrid | 2,903,183 | 70.4% | 69.7% | 93.7% | 96.6% | 92.9% |
| merged only | 2,823,194 | 68.5% | 69.1% | 93.6% | 98.0% | 93.6% |
| 2x150 | 2,985,412 | 73.7% | 65.9% | 93.7% | 95.7% | 93.8% |

With v0.1.0 as released, scored uncapped, the gap rule named 84.1 /
52.4 / 62.4 / 37.3% (2x250 merged, 2x150, for each set, on the earlier
primer trimming); the changes are below. Before exact matches were
named at 2 edits clear (§6), it named 90.4 / 91.2 / 87.0% (Francis &
Šlapeta hybrid, merged, 2x150) and 87.7 / 88.2 / 87.9% (Gilleard).

- **Primer trimming.** In half the Gilleard libraries read 2 starts
  inside the 28S stretch, without NC2 (49% of read pairs). A cutadapt
  trim requiring NC2 discarded them, so the earlier figures for that set
  rested on 41 of 82 libraries, and its four chimera-heavy samples were
  among them. `bin/trim_primers.py` finds NC1 on read 1 (<= 3
  mismatches) and cuts read 2 at the conserved 28S motif
  ATGCTTAAGTTCAGCGGG, 7 bases past NC2 when it is present: every read 2
  starts on the same base either way. 93.9% of Gilleard and 97.9% of
  Francis & Šlapeta pairs pass; chimera removal now drops 0.6% of
  Gilleard reads, from 1.7%. `bin/trim_28s.py`'s anchor stops at the
  motif's end for this: carrying the 4 bases past it, its fit slid into
  ITS2 on reads cut there and cost each 2x150 window 4 edits.
- **Join 2x250 pairs hybrid.** `mergePairs` needs 12 bases of overlap,
  so an ITS2 over about 400 bases never merges at 2x250 and the pair is
  dropped. *Metastrongylus apri* is 413-422 bases: two of the four
  Gilleard pig libraries kept 6 and 665 reads merged and 22,932 and
  19,801 hybrid, and *M. apri* rises from 1.6% to 3.0% of the set's
  reads. Hybrid also keeps the 4.6% of Francis & Šlapeta pairs that
  overlap but fail merging on mismatches. It costs the gap rule 0.7-0.9
  points. Concatenating every pair costs it more where the ITS2 is short:
  *Oesophagostomum radiatum* ITS2 (~221 bases) is read through by read 1
  into 28S, and the window scoring (§6) charges the overrun, 5-6 edits,
  so its reads stop at genus and the gap rule names 68.0% of Gilleard
  reads, against 92.9% hybrid. The k-mer classifiers on cut ASVs are
  unaffected (98%). Cutting the 28S overrun off read 1's window, as
  read 2's is cut, would remove this (§8).
- **Against the published analysis.** Francis & Šlapeta classified with
  assignTaxonomy against their own local database (>96% identity, all 53
  ASVs to species). Pooled by host, every species they report in cattle
  and sheep is within 1-2 points of ours with IdTaxa or assignTaxonomy
  on cut ASVs (cattle *H. placei* 43.5% theirs, 42.9% ours; *O.
  ostertagi* 26.2, 25.7; sheep *H. contortus* 87.5, 88.1). Goats differ
  by 5-7 points (*T. circumcincta* 50.9, 45.3; *Oe. venulosum* 25.0,
  32.1), perhaps from pooling per sample rather than per read. Their two
  *T. circumcincta* lines separate here by length: the goat line, with
  its 6-base insertion, 295 bases and in goats only; the sheep-and-goat
  line 289 bases, mostly in sheep.
- **No arm contradicts another.** Wherever two arms both name a species,
  it is the same species, in all eight runs. These are samples with
  few species, nearly all in the database; the arms differ only in how
  deep they go, and 98-100% of reads reach genus in every arm in the
  table.
- **This database helps both DADA2 classifiers**, mostly through
  curation: with the backbone, *Haemonchus* reads stop at genus because
  of the shared *H. contortus* / *H. placei* variant (§5.3).
- **The 28S stretch must be cut for k-mer classifiers.** ASVs carry
  ~47 bases of 28S that this database's records lack. The span scorer
  ignores them; assignTaxonomy and IdTaxa do not. Uncut, assignTaxonomy
  named 45% of Francis & Šlapeta merged 2x250 reads to species with
  this database; cut, 99.3% (`bin/trim_28s.py`: at
  the 28S anchor, else the last 47 bases). The columns for this
  database above are on cut ASVs.
- **The gap rule names fewer species here** than the k-mer classifiers
  with this database. On v0.1.0, a few weak records decided most of the
  difference: partial records resting on one accession.
  - At 2x250, most of it was *H. placei* reads with an *H. contortus*
    record 3 edits further off (25% of Gilleard reads, 7% of
    Šlapeta's). Typical records of the two species are 4-7 edits apart.
    The 3-edit records were partial copies 1 edit from a full-length
    *H. contortus* record (issue #2). They are now folded into that
    record (§8, item 9), and MK087176, which had been cut on the wrong
    strand, is turned (§3.2).
  - At 2x150, a single partial *H. bedfordi* record
    (`ITS2_2c0973e14d2e`, partial_90) sat 2-3 edits behind both through
    the windows and stopped 49% of Gilleard reads and 37% of Šlapeta's
    at *Haemonchus*. The window scoring gave up to 10 of its first bases
    free; it now charges them (§6), which puts the record +5 or more
    off, as amphid scores it.

  The other reads it stops at genus (6-8% of each set) are mostly an
  exact match to one species with a sister 1 edit off (*Cooperia
  oncophora* / *C. surnabada*; *C. punctata* / *C. spatulata*;
  *Trichostrongylus colubriformis* / *T. capricola*), or < 99% to the
  nearest species with nothing else close (*Oesophagostomum radiatum*
  98.6%, next +29; the goat line of *T. circumcincta* 98.0%, next +7;
  *Oe. venulosum*), or kangaroo *Cloacina*, a genus of many close
  species. Naming exact matches at 2 edits clear (§6) moved 2-6 points
  to species: the rabbit reads (*T. retortaeformis*, *T. vitrinus* +3)
  and the sheep-and-goat line of *Teladorsagia circumcincta*
  (*Ostertagia ostertagi* +2-3; at 2x150 0.5% of Šlapeta reads stay at
  +1). The leave-one-out says the first kind is where assignTaxonomy's
  present-species errors come from, and that the second kind cannot be
  told from a missing relative. The `nearest` and `next` columns put
  both on the record. Two sister pairs look like curation questions
  rather than biology: *C. surnabada* is widely held to be a morph of
  *C. oncophora*, and an *O. ostertagi* record 2-3 edits from *T.
  circumcincta* is closer than two genera should be (§8).

IdTaxa's calls are random draws (bootstraps, and sampling in
training). In the cattle set the large *H. placei* ASVs sit at
confidence 62-76 against a threshold of 60, and across five seeds ASVs
holding 10-23% of reads switched between species and genus; the IdTaxa
columns above are single runs, and the 2x150 / 2x250 differences, and
the backbone IdTaxa column's change since the earlier runs (Francis &
Šlapeta 95.4% then, 81.2% merged now; not tested beyond one seed) are
probably mostly this.

**Which to use.** For samples whose species are all in the database
(livestock), assignTaxonomy or IdTaxa on 28S-cut ASVs with this
database names most, and here agreed with every other arm. Where species may be
missing (wildlife, new hosts, surveys), the gap rule's 3-7x lower
over-classification matters more than its lower recall. amphid should
follow the gap rule (feedback note §5) and report `nearest` and `next`.

### 6.3 Simulated mocks

The simulated nemabiome mocks of Jesudoss Chelladurai et al. (2025,
Parasites & Vectors, 10.1186/s13071-025-07184-1; Zenodo
10.5281/zenodo.17019940): canine (4 species), equine (22) and ruminant
(16), nine samples each, ART-Illumina 2x250 reads simulated from
nemabiome v1.0.0 records, each read named for its species. Run through
`bin/dada2_nemabiome.R` (no primer step: the reads carry none), then
`bin/mock_arms.sh`; scored by `bin/score_mocks.py`
(`results/mocks/`). An ASV's truth is the species of the source record
it copies; ASVs 1-8 edits from a source (simulated errors DADA2 kept)
are scored apart as "noisy".

Every mock species is in the database, so the mocks measure recall,
not over-classification. Exact source copies, read-weighted, species
recall / genus recall. Precision is 100% at species and genus for every
arm in every set: no arm named a wrong species or genus.

| Arm | Canine | Equine | Ruminant |
|---|---|---|---|
| assignTaxonomy x nemabiome | 95.5 / 100 | 100 / 100 | 93.8 / 100 |
| IdTaxa x nemabiome | 100 / 100 | 100 / 100 | 94.6 / 98.2 |
| assignTaxonomy x this db | 95.5 / 100 | 100 / 100 | 85.3 / 99.2 |
| IdTaxa x this db | 100 / 100 | 100 / 100 | 92.8 / 97.7 |
| gap rule x this db | 100 / 100 | 99.0 / 99.0 | 78.9 / 97.8 |

Naming exact matches at 2 edits clear (§6) raised the gap rule from
82.6 / 99.0 / 43.2%: exact copies whose sister is 2-3 edits off
(*C. oncophora*, *C. surnabada* +2-3, 24% of ruminant reads;
*Ancylostoma braziliense* and *A. caninum*, *A. ceylanicum* +2-3,
17% of canine reads). The k-mer rows moved
by 0.1-0.2 points when `bin/trim_28s.py`'s anchor changed (§6.2).
v0.1.0 gave the gap rule 82.6 / 100 / 45.1% species recall. The
ruminant's 2 points are one ASV (1.9% of reads), *T. colubriformis*:
three single-accession partials 1 edit from it were folded into
records 2 edits off (§8, item 9), so *T. retortaeformis* went from
+4 to +3, below the gap. Folding within 1 edit allows exactly this,
and gained 0.7% of reads elsewhere (*Haemonchus*). The equine's 1
point is a mock artefact (below).

Noisy ASVs: species precision 100% in equine and ruminant; in canine
every k-mer arm named 5.7% of them as the wrong species (an *A. caninum*
variant called *A. ceylanicum*: precision 94.3%), and the gap rule
stopped them at genus.

What the mocks can and cannot say:

- **Present species only.** Every source record is in both databases,
  so the mocks re-test §6.1's present case. They reward naming boldly
  and cannot show the gap rule's lower error on missing species. Its
  genus stops are sisters 1 edit apart: *Cooperia oncophora* /
  *C. surnabada* (6.4% of ruminant reads), *Trichostrongylus
  colubriformis* / *T. capricola*, *C. punctata* / *C. spatulata*.
- **They favour the backbone.** The sources are nemabiome's own
  references, boundaries and all; this database's records are re-cut at
  the anchors (§3.2), so ASVs carry flank bases its records lack.
- **The truth carries the backbone's label errors.** The "*Oesophagostomum
  radiatum*" source is MG685662, dropped here as *Mecistocirrus
  digitatus* (§5.2): all five arms call it *M. digitatus*, with
  nemabiome's own database too - hence the paper's "*O. radiatum* not
  detected". Six "*Haemonchus contortus*" sources are the cattle records
  dropped in §5.3; assignTaxonomy (both databases) and the gap rule call
  them *H. placei*.
- **Synonyms count as species.** The mock lists *Cylicostephanus
  poculatus* and *Petrovinema poculatum* as two species, and this
  database holds *Cylicodontophorus* and *Parapoteriostomum mettami*
  apart, which split the k-mer arms' votes on v0.1.0 (4.9% of equine
  reads stopped above genus; none in v0.1.1). Scored as one
  species each.
- **Two sources are not real amplicons.** One *H. placei* source is
  AM410070 as deposited, reverse-complemented, so its ASV (0.75% of
  ruminant reads) is reversed; the classifier now turns it (§6). One
  *Cylicostephanus longibursatus* source is MG738692, a partial record
  ~50 bases short of the ITS2's 3' end. v0.1.0 held it as a record of
  its own and matched it exactly; v0.1.1 folds it into a
  full-length record (§8, item 9), against which the truncated ASV is
  charged for the missing end. That is the equine gap rule's 1.0%
  unassigned. A real amplicon carries the full end.
- **DADA2 at nemabiome.ca's defaults handles ART's reads badly.** With
  error-heavy 3' tails, no truncation and no mismatch allowed in the
  overlap, 3-57% of pairs merged, unevenly: one equine sample kept 4 of
  22 species. Composition error against the simulated truth (15-74%) is
  therefore DADA2's, the same in every arm.

Evaluation stops here. What limits species calls now is curation
(synonyms, shared ITS2 types, missing species), not thresholds or
classifiers.

## 7. Release format

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
- **k-mer classifier files** (`bin/export_kmer.sh`), for the DADA2
  workflow: an assignTaxonomy reference (`bin/export_dada2.py`; a record
  is written only to the rank it reaches), an IdTaxa `trainingSet`
  (`bin/train_idtaxa.R`), and `R/trim_28s.R`, which users must run on
  their ASVs first (§6.2). The records are not padded with a borrowed
  28S stretch to spare that step: most records lack their own, it varies
  by 1-6 bases, and borrowing would invent sequence. `trim_28s.R` gives
  the same cuts as `bin/trim_28s.py` on all 1,143 test and mock ASVs, and
  with the released files reproduced the assignTaxonomy results of §6.2
  exactly.

amphid needs changes before it can use this database:
[feedback/its2_notes_for_amphid.md](feedback/its2_notes_for_amphid.md).

## 8. Open questions

1. **Ties within a family** (27 records stop at family, read whole):
   check as for §5.2? Cross-genus ties are more often real than
   cross-family ones.
2. **Leave-one-out**: done (§6.1, §6.2).
3. **Intraspecific variation**: within one worm, variants differ by
   1-12 edits, as much as some species do (§5.4). The long-read layer
   (§5.5) adds the variants found so far. More long-read assemblies and
   HiFi runs would add more; *H. placei* has none yet.
4. **Contained records of the same species**: done. A partial record
   identical over its span to a full-length record of the same species is
   folded into it (`fold_contained` in `bin/build.py`): 1,115 records in
   v0.1.0, 1,122 in v0.1.1. Its source records stay in `members.tsv`, noted.
   Partials inside another species' full-length record are kept; they
   are real ties.
5. **Shared ITS2 types** (§5.3, §5.4): *Cylicocyclus insigne* /
   *C. ashworthi*; the four look-alike pairs NCBI keeps apart (§5.3); the
   unclear neighbour pairs. Names now follow NCBI (§5.6); refetch
   `curation/ncbi_names.tsv` before each release. These, not thresholds, are what limits species calls
   on missing-species errors (§6.1).
6. **Mock communities**: done (§6.3). They test present species only; a
   missing-species version (each mock species removed in turn) was not
   run.
7. **Unmeasured records** (1,713 in v0.1.1): a lower-identity relative, or an HMM
   for the 28S stretch, could anchor more of them.
8. **ITS1-5.8S-ITS2** for NC5/NC2 nanopore runs, should those be used
   later: nemabiome has a separate database (v0.0.1).
9. **Weak records** (issue #2): done. A weak record - partial or
   unmeasured, resting on one accession - within 1 edit of a robust
   record of its own species is folded into it (`fold_redundant` in
   `bin/build.py`, every pair compared): 752 records in
   v0.1.1, each noted in `members.tsv`. In v0.1.0 a few of them, 2-3
   edits from another species, stopped most of the *Haemonchus* reads
   in §6.2 at genus. A read matching a folded record is at most 1 edit
   further from its host, so a margin can fall by one: in the ruminant
   mock, one *T. colubriformis* ASV fell from +4 to +3 (§6.3).
10. **A record's first bases through 2x150 windows**: done. v0.1.0's
    window scoring left up to `SLACK` (10) of a record's outermost bases
    uncounted; amphid charged them. Both now charge them, and the second
    window is cut at the 28S anchor rather than at a fixed 47 bases
    (§6). It costs 2 points of 2x150 leave-one-out recall and saves a
    quarter of the over-classification (§6.1).
11. **Margins in v0.1.0's tables** (`build.py` `margins()`) came from the
    top 30 other-species candidates by matching bases. Scored uncapped,
    744 records' whole-read margins are lower, 128 fall below 4 edits and
    3 become ties. Done: the cap is lifted, and v0.1.1 carries the
    corrected margins and resolution tables.
12. **Records cut on the wrong strand**: done (§3.2). 36 records turned.
    LN626259 (*Dirofilaria*) maps equally well both ways and is left as
    deposited; AM410070 keeps a 7-base NC1 remnant with an indel.
13. **Read 1 running into 28S.** Where a 2x250 pair does not merge and
    the ITS2 is short, read 1 reads through it into 28S, and the window
    scoring charges the overrun (§6.2: *Oesophagostomum radiatum*, 5-6
    edits). Hybrid joining avoids it, since such pairs merge; cutting
    read 1's window at the 28S anchor, as read 2's is, would make
    concatenated 2x250 safe too.
14. **Two sister pairs from the real data** (§6.2). *Cooperia
    surnabada* sits 1-2 edits from *C. oncophora* and is widely held to
    be its morph; merging them, as the synonym pairs of §5.6, could free
    the 6.4% of the ruminant mock's reads and the real *Cooperia* reads
    still at genus (1 edit off). An *Ostertagia ostertagi* record 2-3
    edits from *Teladorsagia circumcincta* held the sheep-and-goat line
    at genus until exact matches were named at 2 edits clear.
    Check both records' sources before treating either as biology.
