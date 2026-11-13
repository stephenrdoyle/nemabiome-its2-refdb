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

### 2.1 GenBank since the backbone (`bin/genbank_update.py`)

The backbone is a snapshot. `bin/genbank_update.py` fetches the nematode
records it lacks, by four queries: ITS2 records published since 17
March 2026 (1,773), 5.8S/28S records published since (1,919), records
of any date annotated with 5.8S and 28S but never ITS2 (128), and
accessions listed by hand in `curation/genbank_named.tsv` (98: strongylid
species and informal lineages found missing by checking NCBI's
Strongyloidea ITS2 records against the build). Of 3,917, 3,152 sit at or
below a formally named species (informal "sp." names are kept only when
listed by hand) and 1,834 hold an ITS2 region (>= 80 bases at >= 70% to
a build record, by minimap2); the best hit, scored as matches less
twice the differences, stands in for the HMM boundary. `merge` appends
them to the backbone's sequences and taxonomy with a `source` column.

`bin/its2.py` cuts them like any record, but the build takes an update
record only where it was cut at an anchor or its own end: 1,057 were
not, nearly all plant parasites whose record runs on into ITS1 past
where the stand-in boundary lands. 733 update records are kept; with the
2026-11-11 curation pass v0.1.4 has 11,288 records and 1,874 named
species (v0.1.3: 11,135 and 1,864); with the genome scans of 2026-11-13
(§5.5), 11,329 and 1,873. The update is screened before a
build (§5.3): by `bin/species_coherence.py`, and for new records within
2 edits of another species.

Some new records crowd a sister species and are kept as real: a complete
*Trichostrongylus retortaeformis* (KC521376, KC521392, PV439847) 4 edits
from *T. axei*, 3 through 2x150 windows; *Heligmosomoides bibullosus*
1 edit from *H. vandegrifti*; *Metastrongylus asymmetricus* 1 from
*M. salmi*. *Dirofilaria asiatica* was 0-1 edits from eight Sri Lankan
"*D. repens*" records (ON332065-71, PV643877); the two species' genomes
show those are *D. asiatica*, and they are relabelled (§5.5).

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
- **inner NC1**: an NC1 (<= 2 edits) more than 25 bases inside a cut is
  the real start: the start was placed early, by the HMM or an NC1-like
  stretch further back, and 5.8S or ITS1 came with it (KR259650,
  *Ostertagia antipini*: 656 bases, nearly all ITS1).
- **28S stub**: a record ending within 60 bases past the boundary, in a
  28S stub too short for the anchor, is cut at the stub's head as if the
  record ended at the boundary (`part`). PZ347035 and PZ347037, complete
  *Oesophagostomum radiatum* ITS2s with 23 bases of 28S, were kept at
  the HMM's boundary before.
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
  share. 291 records read whole, 429 through 2x150 windows (v0.1.4;
297 and 416 in v0.1.2; 343
  and 473 in v0.1.1, 410 and 522 in
  v0.1.0, 732 and 745 before curation and folding).
- **`species_status`** (`taxonomy.tsv`): `named` if the record reaches
  species; `shared` if it stops above species because its members name
  two or more species (80 records, 95 in v0.1.1: real shared ITS2s, such as *Cylicocyclus
  insigne* / *C. ashworthi*);
  `unnamed` if no member names a species at all (51 in v0.1.4, 32 in
  v0.1.3: hybrids, cf., environmental samples, and the informal "sp."
  names of the GenBank update, §2.1, whose name is in `lineage`).
  "sp."/"spp." count as no species name, as cf. and hybrids do: the
  backbone has none, but 19 update records do, and named as species each
  would have counted as a species of its own in the gap. The classifier's gap counts a shared record as
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

**Lineages.** *Teladorsagia circumcincta* has two ITS2 lines: a
sheep-and-goat line (61 records, mostly 249 bases) and a goat line
(255 bases; a 6-base insertion and 4 substitutions, 10 edits from the
other). Francis & Šlapeta (2023) found the goat line in goats only on
a pasture shared with sheep. Two records here hold it: AJ577464 and
MG651890 (one ITS2; genotype Mdh2C of Leignel et al., and an aoudad in
Tunisia) and AF351782. Rather than split the species,
`curation/lineages.tsv` gives these records `lineage` = "goat line"
and the classifier reports it (§7). Their Australian goat-line ASVs
are 1-5 edits from these, so 8.6% of goat reads stop at genus under the
gap rule (98.0-98.8%) but still report the line.

**The 2026-11-11 curation pass** (`docs/curation_pass_2026-11-11.md`,
53 rows added to `curation/overrides.tsv`) went after the two causes of
the leave-one-out's missing-species errors that a threshold cannot
reach. *One species under two names*: 27 relabels, mostly minor male
morphs of polymorphic species (*Cooperia surnabada* and *C. spatulata*
into *C. oncophora* and *C. punctata*, *Marshallagia occidentalis* into
*M. marshalli*, *Teladorsagia davtiani* into *T. circumcincta*,
*Ostertagia arctica* into *O. gruehneri*). *One odd record*: 8 drops of
records that sit next to another species while the rest of their own
species sits far off (`bin/outlier_records.py`), and 18 drops of
shared ITS2s whose sequence is typical of one of their two names and
not the other (`bin/shared_triage.py`). *Ascaris lumbricoides* and
*A. suum* stay apart by decision. The build has 11,135 records and
1,864 named species (11,146 and 1,869 before); leave-one-out
over-classification fell from 6.4% to 5.9% and species recall rose from
53.4% to 55.3% (§6.1); the ruminant mock gained 9 points, all
*Cooperia* (§6.3); real data gained 1-2 points, also *Cooperia* (§6.2).
Records to check before deciding are listed in the pass's notes.

**After v0.1.3** (the pass's notes, F-H): `bin/species_coherence.py`
lists every named record far from its own species (`own`, edits to the
nearest record of its species) and near another (`other`): *displaced*
(other <= 3, own >= 10), *stray* (own >= 20 and >= 20 from the species'
best-supported record, nothing within 10), *closer*
(`results/curation/species_coherence.tsv`). 65 more drops and 2 relabels
followed, from it, from checking marginal species against every NCBI
record, and from screening the GenBank update (§2.1). Among them:
*Ashworthius sidemi* OP787564 holds an *H. contortus* ITS2; 24 West
African "*H. contortus*" (23 from goats and sheep in Ghana, PZ4068xx,
and PZ569106 from Nigeria) hold the *H. placei* ITS2 that no *H.
contortus* genome copy carries (PZ569106 alone held the main *H.
placei* ASVs of PRJNA1524548 at genus, §6.2); "*Cloacina* sp. 2
AJ-2019" MK573551 is identical to *Zoniolaimus occidentalis*;
"*Trichostrongylus* sp. RR-2017" MG711656 to *Cooperia pectinata*;
*Spiculopteragia asymmetrica* OR284985 and AJ251125 are relabelled
*S. quadrispiculata*: *S. asymmetrica* is its minor male morph. Two pairs of informal
*Oesophagostomum* isolates with identical ITS2s share one `lineage`
name. Dropping partial records (< 75% of the ITS2) of species that have
a full-length record was tested (`--drop-partials`, 394 records, no
species lost) and not adopted: leave-one-out species recall rose 0.3
points and nothing else moved, and 277 of those partials are >= 4 edits
from their species' full-length records - variants a read can match,
which a test of full-length queries cannot count.

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

**The 2026-11-13 scans.** Sixteen more assemblies (ENA and WormBase
ParaSite genomes not scanned before; 47 assemblies in all), chosen
where a species was thin, had a sister problem, or had none: 74 copies
of 15 species. 25 match a backbone record, all of their own species; 49
are new records of 12 species (*D. immitis* 18, *Acanthocheilonema
viteae* 9, *Dictyocaulus viviparus* 6, *Trichuris trichiura* and
*Thelazia callipaeda* 3 each). No species enters. Three settle curation questions:

- ***Dirofilaria repens* / *D. asiatica*.** The *D. repens* genome
  (GCA_008729115; GCA_966211175 is its *Wolbachia*) and the
  *D. asiatica* genome (GCA_053573785) are 14-21 edits apart. The
  Sri Lankan "*D. repens*" records are 1-3 edits from the
  *D. asiatica* genome and are relabelled.
- ***Trichostrongylus capricola*.** JF276023, the only record of the
  species (goat, Iran), is 1 edit from 692 of the 766 ITS2 copies of the
  *T. colubriformis* genome GCA_963978875: inside that species' own
  spread. Dropped; *T. capricola* leaves the database, and the
  *T. colubriformis* reads it held at genus are named (§6.2, §6.3).
- ***Aelurostrongylus abstrusus* GCA_060414415** is *Troglostrongylus
  brevior*: 0-5 edits from *T. brevior* records, 81-192 from
  *A. abstrusus*. Not used.

Twenty assemblies gave no ITS2: all nine *Trichinella* (most likely a
limit of the scanner in Trichinellida) and eleven
short-read assemblies (*Haemonchus placei*, *Trichuris suis*,
*Baylisascaris ailuri*, *Toxascaris*, *Enterobius*, *Syphacia*,
four *Strongyloides*, *Anisakis simplex*).

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
- **Repeat-aware** (since v0.1.3): a difference in how many units a
  simple repeat has - a homopolymer of >= 4 bases, or >= 3 copies of a
  2-3 base unit - costs one edit, not one per base. Each reference is
  scored plain and with its repeats collapsed (two copies of the unit and
  a marker for the copy number), and the better identity is kept. Only
  the best match is scored so: the nearest other species is measured
  plain, so sister species that differ mainly in a repeat's length keep
  their margin. A repeat cut short where a read window stops is not
  charged. Goat *Oe. venulosum* with an (AT)n 16 bases longer than any
  record cost 13-21 edits and held 2.2% of a goat set's reads at family
  (§6.2). Collapsing both sides instead dropped *T. trichiura*,
  *O. volvulus*, two *Baylisascaris* and *Pulchrascaris australis* to
  genus in leave-one-out; collapsing without keeping the plain score
  made a substitution inside a repeat cost several edits. `ITS2_REPEATS=0`
  in the environment turns it off.
- **The gap rule**: species when >= 99% identical and the nearest
  reference of any other species is >= 4 edits further off, or when the
  match is exact and it is >= 2 edits off (sister species share
  variants: an exact match with another species 1 edit off is right
  60% of the time, §6.1); genus at >= 97%, family at >= 95%. `shared` records count as another species,
  `unnamed` ones do not (§5).
- **The tiered rule** (since v0.1.3): the gap rule's species calls, and
  also species at >= 98% with the same 4-edit gap, flagged `tentative`.
  It names intraspecific variants the 99% floor leaves at genus (the
  goat line of *T. circumcincta*, §5.3) at the cost of more names for
  missing species (§6.1). `gap` stays the confident call.
- **What it reports**: identity, ties and gap, and two columns that say
  how close a call that stopped short came. `nearest` is the tied
  species and their identity; `next` is the nearest other species and
  how many edits further it is. A shared record is shown with the
  species its members name (*Cylicocyclus* shared: ashworthi / insigne).

### 6.1 Leave-one-out

As in vespa-18s-refdb. `bin/its2_leave_one_out.py`: up to 5
full-length species-level records per species, 3,230 queries from 1,236
species (after the 2026-11-13 genome scans, §5.5; 3,224 in v0.1.4, with
the curation pass and the GenBank update, §2.1; 3,182 from 1,220 in
v0.1.3, 3,195 from 1,224 before the pass), read whole and through 2x150 windows, as perfect reads.

- **Missing**: every record holding a member of the query's species is
  removed. Any species named is wrong; genus or above is right.
- **Present**: only the query's own record is removed, for species with
  another record.

Reported as recall (right calls / queries), precision (right calls /
calls made at that rank or deeper; a species call is right at genus if
its genus is) and over-classification (a species named when the species
is missing; at genus, a genus named when the genus is missing too).
Specificity is not used: with 1,220 species every query is a negative
for all but one, and every method scores above 99.9%.
(`bin/loo_metrics.py`, `results/leave_one_out/metrics.txt`.)

| Rule (whole) | Species recall | Species precision | Genus recall | Genus precision | Over-classification, species | Over-classification, genus |
|---|---:|---:|---:|---:|---:|---:|
| amphid's V4 limits (99.7% / 99.5% / 95%) | 19.3% | 95.0% | 51.5% | 99.3% | 1.8% | 1.1% |
| exact: species only at 0 edits | 11.5% | 94.6% | 54.3% | 99.3% | 1.2% | 1.1% |
| gap (>= 99% and 4 edits clear, or exact and 2 clear) | 57.3% | 99.3% | 89.4% | 99.4% | 5.4% | 4.5% |
| tiered (gap, plus tentative at >= 98% and 4 clear) | 68.2% | 99.0% | 89.4% | 99.4% | 9.0% | 4.5% |

All repeat-aware (§6). At 2x150 the gap rule gives 60.1%, 99.5%, 90.7%,
99.4%, 5.9% and 4.5%, the tiered rule 68.8%, 99.2% and 9.4% at species
(`results/leave_one_out/`). v0.1.4 gave the gap rule 57.3%, 99.2% and
5.5% whole, the tiered rule 68.2%, 98.9% and 9.3%: the genome scans
(§5.5) cost no recall and remove the *T. capricola* / *T. colubriformis*
errors and most of *D. repens* named as *D. asiatica*. Before repeat-aware scoring the gap rule
gave 56.4%, 99.2% and 5.4% whole, 59.3%, 99.4% and 5.9% at 2x150: it
names 19 more species when present (*Chabertia erschowi*,
*Metastrongylus salmi* and *M. pudendotectus*, *Gongylonema pulchrum*,
*Spirocerca lupi*, *Eucoleus aerophilus* and nematodes of plants and
insects) and loses none; when missing, 10 queries more are named as a
sister (*Baylisascaris*, *Dirofilaria repens* as *D. asiatica*,
*Pulchrascaris*, *Trichuris*), a few fewer otherwise. v0.1.3 gave 55.3%, 98.8% and 5.9% whole,
58.1%, 99.0% and 6.5% at 2x150; the curation pass (§5.3) had moved these from
53.4%, 98.5% and 6.4% whole, 56.3%, 98.7% and 7.0% at 2x150. Before exact matches were named at 2 edits
clear, it gave 52.3%, 98.6% and 6.3% whole, 54.0%, 99.0% and 6.9% at
2x150 (`bin/loo_gap_split.py`, `results/leave_one_out/gap_split.txt`). v0.1.0, scored uncapped, gave 52.6% recall
and 6.5% over-classification read whole, 56.1% and 9.0% at 2x150: the
stricter window scoring costs 2 points of recall at 2x150 and saves a
quarter of the over-classification. The V4 and exact rules seldom
name a species even when it is there. The gap rule names three to five
times as many, and is wrong less often on present species, because the
margin to the next species is what separates right from wrong, not the
identity. Its missing-species errors (175 of 3,230 queries read whole;
178 in v0.1.4, 173 before repeat-aware scoring, 187 in v0.1.3, 204 before the curation pass) are of three kinds. Most (142) are 99-100% to a sister species with
nothing else within 4 edits: two species under 1% apart, one of them
missing. 33 are identical to a record of another species, a shared ITS2
or a mislabelled record. 22 of the 175 land in another genus, nearly
all curation questions (§5.3). No threshold catches any of them;
curating names and records can (`docs/curation_pass_2026-11-11.md`).

A third of present-species reads stop at genus. 27% of queries are
below 99% identity to their own species (intraspecific variation; §5.4)
and 18% have a gap of 0-3 edits (sister species whose ITS2 barely
differs). Could a rule per genus or family do better?
`bin/loo_group_rules.py` let each genus or family choose its own
identity floor (97-99%) and gap (1-10 edits), chosen on four fifths of
the species and scored on the rest. On held-out species, per-group rules
recovered 2.0-3.6 points of right calls for 0.3-0.5 more present-wrong
and 1.2-1.8 more missing-wrong. Loosening the global rule trades 11.3
points for 0.8 and 3.6 (`results/leave_one_out/group_rules.txt`). No
per-group rule beats the global one by more than fractions of a point,
from five splits of 1,220 species (v0.1.3; the rerun on `main` gives the
same picture); one rule is kept.

The identity floor is the larger part of what the gap rule gives up.
Of the present-species queries it stops at genus, half have the next
species >= 4 edits off but are under 99% to their own; assignTaxonomy
names 91% of these right. Lowering the floor with the gap kept at 4
(whole, before repeat-aware scoring):

| Floor | Present: right | Missing: wrong |
|---|---:|---:|
| 99% (gap rule) | 56.4% | 5.4% |
| 98.5% | 63.4% | 7.2% |
| 98% | 67.7% | 9.0% |
| 97% | 70.5% | 14.1% |
| nearest species >= 1 edit ahead, >= 97% | 80.2% | 21.7% |

IdTaxa (80.3%, 20.5%) behaves like the last row; assignTaxonomy (86.5%,
42.6%) names a species almost whenever one is near. 98% is the tiered
rule's floor: on the mocks it names no wrong species, and on real data
it names the reads the gap rule left at genus that look like
intraspecific variation (§6.2). Its calls are flagged rather than
merged into the gap rule's because a query 1-2% from its own species
looks the same as one 1-2% from a missing sister.

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
| assignTaxonomy | 86.5% | 95.6% | 97.7% | 99.4% | 42.6% | 24.5% |
| IdTaxa | 80.3% | 97.2% | 94.0% | 99.5% | 20.5% | 9.1% |
| gap rule | 57.3% | 99.3% | 89.4% | 99.4% | 5.4% | 4.5% |

At 2x150: assignTaxonomy 87.3% recall, 96.2% precision, 38.0%
over-classification; IdTaxa 80.6%, 97.8%, 20.1%; the gap rule 60.1%,
99.5%, 5.9%. The three sit on one trade-off line. At genus they are
equally precise; at species they differ in how often they name a
species that is not there. assignTaxonomy never declines: a read from a
missing species gets a species name 44% of the time. IdTaxa halves
that. The gap rule over-classifies least, by 3-7 times, and names
fewest. If 20% of a sample's reads come from species missing from the
database, the share of species calls that are right is about 85%
(assignTaxonomy), 90% (IdTaxa) and 96% (gap rule). Of missing-species
reads whose genus is still in the database, the gap rule names the
genus for 26% (precision 94%), IdTaxa 40% (94%), assignTaxonomy 71%
(94%); most of the gap rule's rest stay unassigned. Reads left
unassigned by IdTaxa and the gap rule are not junk: they are reads whose
nearest reference is too far for a family (< 95%).

**Where the k-mer classifiers' extra species calls come from**
(`bin/loo_kmer_vs_gap.py`, `results/leave_one_out/kmer_vs_gap.txt`;
whole reads, points of all present queries, right / wrong species calls
on the reads the gap rule stops short on):

| Why the gap rule stops short | Share | assignTaxonomy | IdTaxa |
|---|---:|---:|---:|
| own species best, but < 99% | 22.5 | 20.5 / 0.3 | 18.5 / 0.1 |
| own species >= 99%, a sister 0-3 edits off | 8.9 | 6.5 / 0.5 | 4.8 / 0.3 |
| another species closest | 8.1 | 2.1 / 2.5 | 0.5 / 1.5 |
| tied with another species | 2.1 | 0.9 / 0.4 | 0.0 / 0.0 |

- **The 99% floor is most of it.** Two thirds of assignTaxonomy's
  29-point recall lead are reads from a divergent variant of a present
  species, which it names almost always correctly. Its confidence is
  relative (does this species keep winning the bootstrap?), the gap
  rule's absolute (is the read within 1% of a record?). The same
  relative judgement names missing species: of hidden queries whose
  nearest record is another species, assignTaxonomy names a species for
  35 points of all hidden queries, IdTaxa for 16.
- **One or two diagnostic sites suffice for them.** The gap rule wants 4
  edits clear of the next species (2 for an exact match); the k-mer classifiers resolve sisters
  1-3 edits apart, mostly correctly when the species is present. The
  real-data *Cooperia* genus stops were this kind before the curation
  pass folded *C. surnabada* into *C. oncophora* (+1).
  Tuned on held-out species, the gap is chosen at 2-4 edits for whole
  reads and at 3 or 4 at 2x150, within a point of the current rule (`results/leave_one_out/scaled_floor.txt`,
  plain rule, lam 4).
- **Where another species is closest, they gain nothing**: assignTaxonomy
  is right as often as wrong, IdTaxa wrong more often than right.

They do not see what identity misses. Within one identity band of the
best hit, the k-mer classifiers name a species for a present species'
variant more often than for a missing species, but not by much: at
98-99%, assignTaxonomy 94% against 67% (333 right calls to 158 wrong),
IdTaxa 90% against 56%; below 98%, wrong calls outnumber right ones.
About 2 right calls per wrong one is also what a distance-scaled floor
gives the gap rule: below 99%, a species at >= 98% when the next species
is >= max(8, 3 x edits to the best) further off. Chosen on four fifths
of the species and scored on the rest, as for the per-group rules
(§6.1), it adds 9 points of recall at unchanged present-species
precision for 2.8 more points of over-classification, and is chosen
only when a wrong call is valued at 2 right ones or fewer; at 4 or 8
every fold keeps the current rule (`bin/loo_scaled_floor.py`,
`results/leave_one_out/scaled_floor.txt`). So the difference is where
each sets its bar, not information: the k-mer classifiers are tuned for
the species being in the database, the gap rule for it perhaps not
being. The current rule is kept, with no relaxed mode.

**Real data**: two NC1/NC2 sets from the Šlapeta group, one published,
sequenced 2x251.
Francis & Šlapeta (2023, Parasitology 150:672-682, 10.1017/S0031182023000380; PRJNA936713): larval cultures
from 10 cattle, 10 feral goats and 10 sheep on one New South Wales
property. PRJNA1524548 (Šlapeta group, University of Sydney; registered 2026, unpublished): 82 libraries, cattle larval cultures
and faecal samples from kangaroos, pigs and rabbits, two libraries per
sample. Each run from raw reads by `bin/real_data.sh`: primer trimming
(`bin/trim_primers.py`), DADA2 denoised once and joined three ways
(merged, concatenated, hybrid), and the same reads cut to 2x150. Share
of reads named to species (`results/dada2_comparison/`, with an HTML
report per set):

| Set | Reads | assignTaxonomy x nemabiome | IdTaxa x nemabiome | assignTaxonomy x this db | IdTaxa x this db | gap rule x this db | tiered x this db |
|---|---:|---:|---:|---:|---:|---:|---:|
| Francis & Šlapeta 2x250, hybrid | 488,062 | 72.0% | 81.7% | 99.5% | 93.2% | 95.7% | 99.4% |
| merged only | 465,576 | 70.2% | 81.2% | 99.5% | 93.0% | 96.2% | 99.4% |
| 2x150 | 484,432 | 76.0% | 88.2% | 99.1% | 99.3% | 95.2% | 98.5% |
| PRJNA1524548 2x250, hybrid | 2,903,183 | 70.4% | 69.7% | 95.1% | 74.8% | 94.3% | 98.3% |
| merged only | 2,823,194 | 68.5% | 69.1% | 99.1% | 75.4% | 94.7% | 98.5% |
| 2x150 | 2,985,412 | 73.7% | 65.9% | 94.8% | 80.3% | 94.7% | 98.3% |

After the 2026-11-13 genome scans (§5.5), repeat-aware (§6).

- **The genome scans** change the gap rule only through *T. colubriformis*,
  named now that *T. capricola* JF276023 is dropped: 2,052 Francis &
  Šlapeta reads (2,012 of them goats; goats 87.2% to 88.6% by the gap
  rule, 97.5% to 98.9% by the tiered rule) and 164 PRJNA1524548 rabbit
  reads. v0.1.4 gave the gap rule 95.3 / 95.8 / 94.8% and the tiered
  rule 99.0 / 99.0 / 98.1% on Francis & Šlapeta (hybrid, merged,
  2x150); PRJNA1524548 is unchanged. No read changed from one species
  to another. IdTaxa x this db moves with its training (below): v0.1.4
  gave 99.8 / 99.9 / 99.7% and 75.9 / 75.2 / 85.5%.

- **Repeat-aware scoring** changes only Francis & Šlapeta at 2x250
  (93.9 / 94.7% hybrid / merged before). Goat *Oe. venulosum* ASVs carry
  an (AT)n 16 bases longer than any record; their pairs do not merge, and
  concatenated they were 13-21 edits from Y10790 and held 2.2% of goat
  reads (0.9% of the set) at family or unassigned. Now 1.3% of the set
  is named *Oe. venulosum* and the rest reach genus; no read is left at
  family in any host. The repeat sits outside both 2x150 windows, so
  2x150 is unchanged. With the full-quality reads (SRA's normalized
  files; ENA and SRA Lite carry Q30/Q3 only) the gap rule names 94.1% of
  hybrid reads before and 95.4% after; quality trimming (cutadapt -q 20,
  25) or fixed truncation (200/180, 180/170) did not help.
- **The tiered rule** names 98-99% of every set but PRJNA1524548
  concatenated (94.6%). The tentative calls are the goat line of
  *T. circumcincta* (98.0-98.8% to the only goat-line records, from
  Europe and Tunisia; 2.4% of Francis & Šlapeta), *Oe. venulosum* (3-4
  edits; 1.4%), cattle *Oe. radiatum* (3 edits; 0.4% and 3% of
  PRJNA1524548), and kangaroo *Cloacina* (98.0-98.7%; 0.5%), the one
  group where they may be missing sisters rather than variants: 72
  *Cloacina* species are named here, of many more described. Goats go
  from 83.1% to 97.5% to species, kangaroos from 62.6% to 76.9%; 15-17% of
  kangaroo reads stay at family or unassigned under every rule.

v0.1.3 gave the gap rule 93.8 / 94.7 /
94.5% (Francis & Šlapeta hybrid, merged, 2x150) and 95.2 / 95.7 / 95.9%
(PRJNA1524548). The update's one regression was a mislabel: PZ569106, a
Nigerian "*H. contortus*" with the *H. placei* ITS2, held the main *H.
placei* ASVs at genus (95.9% to 71.2% of PRJNA1524548 2x150 reads) until
dropped (§5.3). What remains is two kinds, both more honest than before:
a complete *Oe. radiatum* record (PZ347037, §3.2) now holds the partials
the old calls rested on, and scored against its full length some ASVs
are 98.3-98.8% (1 point of PRJNA1524548 to genus); and a complete
*T. retortaeformis* record 4 edits from *T. axei* is 3 off through 2x150
windows (§2.1). The IdTaxa x this db column on PRJNA1524548 is one
unlucky training: seed 1 leaves the main *H. placei* ASVs at genus;
seeds 2 and 3 name 99.0% and 99.2% of hybrid reads (below).

With v0.1.0 as released, scored uncapped, the gap rule named 84.1 /
52.4 / 62.4 / 37.3% (2x250 merged, 2x150, for each set, on the earlier
primer trimming); the changes are below. Before exact matches were
named at 2 edits clear (§6), it named 90.4 / 91.2 / 87.0% (Francis &
Šlapeta hybrid, merged, 2x150) and 87.7 / 88.2 / 87.9% (PRJNA1524548); before the
curation pass (§5.3), 92.0 / 92.9 / 92.4% and 92.9 / 93.6 / 93.8%. The
curation gain is *Cooperia*, named now that its male morphs are folded
into *C. oncophora* and *C. punctata* (§5.3).

- **Primer trimming.** In half the PRJNA1524548 libraries read 2 starts
  inside the 28S stretch, without NC2 (49% of read pairs). A cutadapt
  trim requiring NC2 discarded them, so the earlier figures for that set
  rested on 41 of 82 libraries, and its four chimera-heavy samples were
  among them. `bin/trim_primers.py` finds NC1 on read 1 (<= 3
  mismatches) and cuts read 2 at the conserved 28S motif
  ATGCTTAAGTTCAGCGGG, 7 bases past NC2 when it is present: every read 2
  starts on the same base either way. 93.9% of PRJNA1524548 and 97.9% of
  Francis & Šlapeta pairs pass; chimera removal now drops 0.6% of
  PRJNA1524548 reads, from 1.7%. `bin/trim_28s.py`'s anchor stops at the
  motif's end for this: carrying the 4 bases past it, its fit slid into
  ITS2 on reads cut there and cost each 2x150 window 4 edits.
- **Join 2x250 pairs hybrid.** `mergePairs` needs 12 bases of overlap,
  so an ITS2 over about 400 bases never merges at 2x250 and the pair is
  dropped. *Metastrongylus apri* is 413-422 bases: two of the four
  PRJNA1524548 pig libraries kept 6 and 665 reads merged and 22,932 and
  19,801 hybrid, and *M. apri* rises from 1.6% to 3.0% of the set's
  reads. Hybrid also keeps the 4.6% of Francis & Šlapeta pairs that
  overlap but fail merging on mismatches. It costs the gap rule 0.7-0.9
  points. Concatenating every pair costs it more where the ITS2 is short:
  *Oesophagostomum radiatum* ITS2 (~221 bases) is read through by read 1
  into 28S, and the window scoring (§6) charged the overrun, 5-6 edits,
  against the partial records that were all this database held: the gap
  rule named 68-69% of PRJNA1524548 reads concatenated, against 92-95%
  hybrid, until v0.1.3. With a complete *Oe. radiatum* record (PZ347037,
  §3.2) the window meets its own record's end and the gap rule names
  92.4%. The k-mer classifiers on cut ASVs were never affected (99%).
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
    record 3 edits further off (25% of PRJNA1524548 reads, 7% of
    Francis & Šlapeta's). Typical records of the two species are 4-7 edits apart.
    The 3-edit records were partial copies 1 edit from a full-length
    *H. contortus* record (issue #2). They are now folded into that
    record (§8, item 9), and MK087176, which had been cut on the wrong
    strand, is turned (§3.2).
  - At 2x150, a single partial *H. bedfordi* record
    (`ITS2_2c0973e14d2e`, partial_90) sat 2-3 edits behind both through
    the windows and stopped 49% of PRJNA1524548 reads and 37% of Francis & Šlapeta's
    at *Haemonchus*. The window scoring gave up to 10 of its first bases
    free; it now charges them (§6), which puts the record +5 or more
    off, as amphid scores it.

  The other reads it stops at genus (4-6% of each set, 6-8% before the
  curation pass) are an exact match to one species with a sister 1 edit
  off (*Trichostrongylus colubriformis* / *T. capricola*, until
  JF276023 was dropped, §5.5; until the
  curation pass folded the male morphs, *Cooperia oncophora* / *C.
  surnabada* and *C. punctata* / *C. spatulata*, §5.3), or < 99% to the
  nearest species with nothing else close (*Oesophagostomum radiatum*
  98.6%, next +29; the goat line of *T. circumcincta* 98.0%, next +7;
  *Oe. venulosum*), or kangaroo *Cloacina*, a genus of many close
  species. Naming exact matches at 2 edits clear (§6) moved 2-6 points
  to species: the rabbit reads (*T. retortaeformis*, *T. vitrinus* +3)
  and the sheep-and-goat line of *Teladorsagia circumcincta*
  (*Ostertagia ostertagi* +2-3; at 2x150 0.5% of Francis & Šlapeta reads stay at
  +1). The leave-one-out says the first kind is where assignTaxonomy's
  present-species errors come from, and that the second kind cannot be
  told from a missing relative. The `nearest` and `next` columns put
  both on the record. One sister pair still looks like a curation
  question rather than biology: an *O. ostertagi* record 2-3 edits from
  *T. circumcincta* is closer than two genera should be (§8). The other,
  *C. surnabada*, a morph of *C. oncophora*, is now folded into it.

IdTaxa's calls are random draws (bootstraps, and sampling in
training). In the cattle set the large *H. placei* ASVs sit at
confidence 62-76 against a threshold of 60, and across five seeds ASVs
holding 10-23% of reads switched between species and genus; the IdTaxa
columns above are single runs, and the 2x150 / 2x250 differences, and
the backbone IdTaxa column's change since the earlier runs (Francis &
Šlapeta 95.4% then, 81.2% merged now; not tested beyond one seed), and
IdTaxa x this db on PRJNA1524548 2x150 (95.7% before the curation pass,
87.2% after, with *H. placei* partly at genus), are probably mostly this.
In v0.1.4 it is shown: trained with seed 1, IdTaxa x this db
names 75.9% of PRJNA1524548 hybrid reads to species, the main *H.
placei* ASVs at genus; trained with seeds 2 and 3 on the same database,
99.0% and 99.2% (v0.1.3, seed 1: 99.1%). After the genome scans the same
happens on Francis & Šlapeta: their main *H. placei* ASV (6.5% of reads)
scored 61.8 in v0.1.4 and now falls just below 60, so IdTaxa x this db
names 93.2% of hybrid reads with seed 1, 93.3% with seed 2 and 99.9%
with seed 3.

**Which to use.** For samples whose species are all in the database
(livestock), assignTaxonomy or IdTaxa on 28S-cut ASVs with this
database names most, and here agreed with every other arm; the tiered
rule names nearly as many (98-99%) with half IdTaxa's
over-classification. Where species may be
missing (wildlife, new hosts, surveys), the gap rule's 3-7x lower
over-classification matters more than its lower recall, and tentative
calls are best read as genus. amphid should
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
| assignTaxonomy x nemabiome | 95.5 / 100 | 100 / 100 | 93.0 / 100 |
| IdTaxa x nemabiome | 100 / 100 | 100 / 100 | 91.5 / 98.3 |
| assignTaxonomy x this db | 95.5 / 100 | 100 / 100 | 88.7 / 99.2 |
| IdTaxa x this db | 98.1 / 100 | 100 / 100 | 95.3 / 99.4 |
| gap rule x this db | 100 / 100 | 99.0 / 99.0 | 96.9 / 99.0 |

After the 2026-11-13 genome scans (§5.5), repeat-aware (§6). v0.1.4
gave the gap rule 90.8 / 99.0 on the ruminant mock: the 6 points are
*T. colubriformis*, named now that *T. capricola* JF276023, 1 edit off,
is dropped. The IdTaxa columns move by a point or two with retraining
(§6.2); v0.1.4 gave 99.7 (canine) and 96.6 / 99.2 (ruminant). Before
repeat-aware scoring the gap rule gave 89.7 / 99.0 on the ruminant
mock. The tiered rule gives the same as the gap rule on the ruminant
mock and, on equine ASVs 1-8 edits off their source, 82.7% against the
gap rule's 73.0%, at 100% precision. v0.1.3 gave assignTaxonomy x this db
88.7 / 99.2, IdTaxa 96.4 / 100 and 95.3 / 96.7, and the gap rule
87.8 / 97.9 on the ruminant mock.

The curation pass (§5.3) raised the gap rule on the ruminant mock from
78.9% and assignTaxonomy x this db from 85.3%: *Cooperia* reads are
named now that the male morphs *C. surnabada* and *C. spatulata* are
folded into *C. oncophora* and *C. punctata*. `bin/score_mocks.py` reads
those whole-species relabels in `curation/overrides.tsv` as synonyms,
for the mock labels and the calls alike, so the 4.3% of ruminant reads
once scored "ambiguous source" (*C. punctata* and *C. spatulata* sources
sharing an ASV) are now clean copies. The nemabiome rows differ from
earlier runs by this scoring alone; the IdTaxa rows move by a few points
between runs (§6.2).

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
  genus stops were mostly sisters 1 edit apart: *Trichostrongylus
  colubriformis* / *T. capricola* (5.8% of ruminant reads) until
  JF276023 was dropped (§5.5); *Cooperia oncophora* / *C. surnabada*
  (6.4%) and *C. punctata* / *C. spatulata* until the curation pass.
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
- `taxonomy.tsv` ends with `lineage`, a published line within a species
  (`curation/lineages.tsv`): blank for nearly every record. The species
  name is unchanged, so the ranks, resolution tables and k-mer files are
  as without it. The classifier reports it in a `lineage` column when
  every tied best record carries it.
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
    concatenated 2x250 safe too. A complete record of the species also
    removes it: with PZ347037 in (§3.2), concatenated PRJNA1524548 rose
    from 69% to 92% of reads named.
14. **Two sister pairs from the real data** (§6.2). *Cooperia
    surnabada*, 1-2 edits from *C. oncophora* and widely held to be its
    morph: done, folded into it by the 2026-11-11 curation pass (§5.3),
    which named the 6.4% of the ruminant mock's reads and the real
    *Cooperia* reads that stopped at genus. An *Ostertagia ostertagi*
    record 2-3 edits from *Teladorsagia circumcincta* held the
    sheep-and-goat line at genus until exact matches were named at 2
    edits clear; check its source before treating it as biology.
15. **The curation pass's open groups**
    (`docs/curation_pass_2026-11-11.md`): records to check before
    deciding (group C; OP787564, the one *Ashworthius sidemi* record,
    is now dropped as an *H. contortus* ITS2, §5.3), the
    *Marshallagia* shared records, and the other `sib=0` records of
    `results/curation/outlier_records.txt` outside livestock and
    companion animals.
16. **After the GenBank update** (§2.1): why *Hsiungia* failed the
    update's ITS2 region test. Settled by the 2026-11-13 genome scans
    (§5.5): *T. capricola* JF276023 lies inside the *T. colubriformis*
    genome's own spread and is dropped; the Sri Lankan "*D. repens*"
    records are *D. asiatica*. Open from the scans: no ITS2 found in any
    *Trichinella* assembly (the scanner, most likely); and
    GCA_060414415, deposited as *Aelurostrongylus abstrusus*, carries a
    *Troglostrongylus brevior* ITS2, which its submitters may want to
    know. Informal lineages are kept only where
    listed by hand: a rule for the rest is open.
17. **Repeat-aware scoring and the tiered rule** (§6): done, as defaults
    beside the gap rule. Open: scoring each reference twice roughly
    doubles window scoring, most felt on concatenated tables; the
    tentative calls on kangaroo *Cloacina* cannot be checked without
    more *Cloacina* references. Australian goat-line *T. circumcincta*
    references (2026-11-12): GenBank has none (the Australian
    *T. circumcincta* ITS2s are 14 sheep fragments, MW161470-83), and
    nothing in it is closer to the Francis & Šlapeta goat-line ASVs
    than AJ577464/MG651890 (1-5 edits). Their own reference (database v2,
    Usyd059-063, and *Oe. venulosum* Usyd065-068) holds the exact
    variants, but as ASVs of the same study named by similarity; they
    are left out, and the goat line stays annotated on our records.
    Australian *Oe. venulosum* references (2026-11-12): GenBank's one
    Australian *Oe. venulosum* record is not ITS2; of 142
    *Oesophagostomum* ITS2s, Y10790 (sheep) is nearest the goat ASVs
    (1-4 edits for the main merged ones), then MG651896 (aoudad,
    Tunisia), both held, and the rest 15 or more edits off. No reference
    carries the expanded repeat: the longest (AT)n is 14 bases in GenBank
    (Y10790, MG651896) and 12-16 in the paper's Usyd065-068, against
    14-38 in the goat reads (about 4,400 reads at 22-38). Repeat-aware
    scoring (§6), not a reference, is what names them.
