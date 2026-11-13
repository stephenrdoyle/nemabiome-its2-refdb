# Curation pass, 2026-11-11

Status: **A, B and E applied** (53 rows added to `curation/overrides.tsv`,
2026-11-11; v0.1.2); C left to check; D: *Ascaris lumbricoides* and *A. suum*
stay apart (decided), the rest open. Rebuilt: 11,135 records, 1,864
named species, 54 source records dropped (11,146, 1,869 and 28 before).
**F, G, H and I applied after v0.1.3** (2026-11-12; released in v0.1.4): 65 more
drops and 2 relabels, the GenBank update, and two fixes to the cut.
Rebuilt: 11,288 records, 1,874 named species, 119 source records
dropped (I, below).

## Effect

Gap rule, before the pass / A and B / A, B and E
(`results/leave_one_out/`, `results/mocks/`, `results/dada2_comparison/`):

| Test | Before | A+B | A+B+E |
|---|---:|---:|---:|
| Leave-one-out, species present, right species (whole) | 53.4% | 54.4% | 55.3% |
| same, 2x150 | 56.3% | 57.2% | 58.1% |
| Leave-one-out, species missing, wrong species (whole) | 6.4% | 6.3% | 5.9% |
| same, 2x150 | 7.0% | 6.7% | 6.5% |
| Ruminant mock, species recall | 78.9% | 87.8% | 87.8% |
| Canine / equine mocks | 100 / 99.0% | same | same |
| Francis & Šlapeta 2x250 hybrid, reads to species | 92.0% | 93.8% | 93.8% |
| PRJNA1524548 2x250 hybrid, reads to species | 92.9% | 95.2% | 95.2% |

A did the mock and real-data work: the *Cooperia* morphs folded, so
*Cooperia* reads are named. E did leave-one-out work only: none of its
18 records sits near a mock or real-data ASV. Missing-species errors
(whole) fell from 204 to 187; 22 still land in another genus. Precision
on present species rose from 98.5% to 98.8%. assignTaxonomy and IdTaxa
on the same queries: species recall 85.4% to 86.0% and 76.9% to 78.8%,
over-classification 43.9% to 43.5% and 22.7% to 23.2% (A+B not run
for them).

## Why

Of the gap rule's 204 missing-species errors in the leave-one-out (read
whole), 161 land on a sister species under 1% away and 43 on a record
identical to the query. 22 of them land in another genus (DESIGN §6.1).
A threshold can't remove any of these. Curation can, where the cause is
one of two things:

- **One species under two names.** Usually a major and minor male morph
  of a polymorphic species. The two names compete, so reads of either
  stop at genus.
- **One odd record.** A record that sits next to another species while
  the rest of its own species sits far away. It holds that species' reads
  at genus, and when the species is missing it gives a wrong name.

`bin/outlier_records.py` finds the second kind from the build's
`species_margins.tsv` (`results/curation/outlier_records.txt`). For each
record it lists the records of another species that have it as their
nearest other species within 3 edits. `sib=0` means no other record of
its own species comes that close. 63 named records have `sib=0`; this
pass covers the livestock and companion-animal ones. Earlier examples of
the same pattern are the *H. contortus* partials folded in v0.1.1 and
the *H. bedfordi* partial (DESIGN §6.2).

## A. Relabel: one species under two names (27 rows)

Each rests on a published study, in most cases the very study that
deposited the record.

| From | To | Records | Evidence |
|---|---|---:|---|
| *Teladorsagia davtiani* | *T. circumcincta* | 1 | Minor male morph (Grillo et al. 2008: no microsatellite differentiation; ITS2 0.8% apart). Treated the same way as *O. trifurcata*, already relabelled. |
| *Marshallagia occidentalis* | *M. marshalli* | 13 | Dallas et al. 2001, *Syst. Parasitol.* 50:101 (AH010211): "conspecific". The records interleave at 0-1 edits. |
| *Ostertagia arctica* | *O. gruehneri* | 1 | Minor morph; ITS1 and ITS2 virtually identical (Dallas et al. 2000). |
| *Cooperia surnabada* | *C. oncophora* | 1 | Newton et al. 1998, *Int. J. Parasitol.* 28:331, the source of AJ000032 itself: "one species". This record is the nearest other species for 56 of 64 *C. oncophora* records. |
| *Cooperia spatulata* | *C. punctata* | 9 | Ramünke et al. 2018, *PLoS ONE*: *C. spatulata* is a morphotype of *C. punctata*. |
| *Heligmosomoides polygyrus* (DQ408624, AY333382) | *H. bakeri* | 2 | Both were deposited as "*H. polygyrus bakeri*", the laboratory lineage. Cable et al. 2006 (DQ408624's own study) show it is a separate species. They sit 0-1 edits from the *H. bakeri* genome copies and 10-14 from wood-mouse *H. polygyrus*. |

Expected effect: the leave-one-out and real data will name *Cooperia*
and *Teladorsagia* reads that now stop at genus. The ruminant mock's
6.4% of *C. oncophora* reads still at genus should reach species.

## B. Drop: one odd record (8 rows)

None of these records has a publication. Each holds another species'
records at genus, while the rest of its own species sits well away from
that species.

| Accession | Deposited as | Problem |
|---|---|---|
| ON497074 | *Ostertagia ostertagi* (sheep/goat, Jammu; partial) | Identical to *T. circumcincta* MN845174. Nearest other species, within 3 edits, for 18 *T. circumcincta* records; the other 21 *O. ostertagi* records sit a median 20 edits from anything. This is the record behind DESIGN §8 item 14. |
| ON004118 | *O. ostertagi* (sheep, Iraq) | Identical to *O. leptospicularis* KC998721 (Bisset et al. 2014). Holds 9 of 19 *O. leptospicularis* records. |
| KT215383 | *Cooperia oncophora* (cattle, Kenya) | Identical to *C. pectinata* MH267780 (Ramünke et al. 2018). Holds all 5 *C. pectinata* records. |
| PX596677, PX596678 | *Toxocara cati* (Iraq; partial) | The isolates are named "TCani3" and "TCani4". They sit 2-3 edits from *T. canis*; other *T. cati* records sit a median 20 away. |
| PX497933 | *Ancylostoma caninum* (dog, Laos; partial) | Within 3 edits of 32 of 38 *Necator americanus* records; none of the other 38 *A. caninum* records is close. |
| MW528471 | *Oesophagostomum bifurcum* (chimpanzee, Cameroon; partial) | Holds 8 of 46 *Oe. stephanostomum* records, which is also a chimpanzee parasite. |
| JQ995321 | *Ascaridia columbae* (dove, Texas) | 1 edit from *A. galli* and holds 9 of 69 *A. galli* records; the other 5 *A. columbae* records sit a median 9 away. |

Relabelling rather than dropping would be defensible for ON497074 (as
*T. circumcincta*) and PX596677/8 (as *T. canis*). Dropping is the
cautious choice: a misidentified record is kept only when the right name
is certain.

## C. To check before deciding (not applied)

- **MK087183** *T. axei* (sheep, Egypt; partial) holds 9 of 10 *T.
  retortaeformis* records. **ON721592** *T. vitrinus* (goat, Bangladesh;
  partial) holds 13 of 41 *T. axei* records. Both look like odd records,
  but 3 and 6 other records of their species are also within 3 edits,
  so *Trichostrongylus* may simply be close. Check alignments.
- **JF276023**, the only *T. capricola* record (goat, Iran), holds 20 of
  48 *T. colubriformis* records. Dropping it removes the species.
  Dropped on 2026-11-13: it is 1 edit from 692 of the 766 ITS2 copies of
  the *T. colubriformis* genome GCA_963978875 (docs/DESIGN.md §5.5).
- **MT193647** *T. colubriformis* (Bactrian camel, China) is identical to
  *T. probolurus* EF427623 (Kuznetsov 2007). *T. probolurus* is a camel
  parasite, so this record may be *T. probolurus*. Dropped in E.
- **OP787564**, an *Ashworthius sidemi* record (European bison;
  partial), holds 7 *H. contortus* records. It is *H. contortus* ITS2:
  proposed in F.
- **PZ103494 / PZ103495** (*Cylicostephanus calicatus* / *Cylicodontophorus
  bicoronatus*, one Malaysian horse study) sit 2 edits apart, while both
  species otherwise sit a median 11 from anything. One of the two is
  probably mislabelled.
- ***Spiculopteragia boehmi*** (AJ577460, MW839588) is identical to *S.
  spiculoptera* (DQ354326, Liénard et al. 2006) and holds 9 of 11 of its
  records. I found no published synonymy; Liénard et al. paired *S.
  spiculoptera* with *S. mathevossiani*, not *S. boehmi*.
- **Cross-genus pairs from the leave-one-out**, already on the DESIGN
  §5.3 list:
  - *Rugopharynx alpha* LN906946 (Chilton et al. 2016) is identical to
    *Papillostrongylus barbatus* KJ776392 (Lott et al. 2014);
  - *Ektaphelenchoides compasi* DQ257622 is 2 edits from *Bursaphelenchus
    thailandae*, from the same lab;
  - *Aspidodera raillieti* / *Nematomystes scapteromi* come from one
    study on host switching (Jiménez et al. 2012), and may be real.

## C2. Found after the rebuild

- Two shared *Marshallagia* records hold *M. marshalli* reads: MN845179
  (*M. marshalli*) is identical to MT110927 (*M. schumakovitschi*), and
  MG011724 / MZ148588 (*M. marshalli*) to MT118028 / MZ148613 (*M.
  trifida*). Kuchboev et al. (MT110967's study) treat *Marshallagia* as
  polymorphic; check whether these are further morphs.

## E. Shared ITS2 that fits one name only (18 rows, applied)

`bin/shared_triage.py` (`results/curation/shared_triage.tsv`) takes every
place in the curated build where one ITS2 carries two species' names:
shared records, and named records 0 edits apart over their span (233
cases). For each name it asks how the sequence sits among that
species' other records: the edits to the nearest one, the share within 3
edits, and the species' own spread (median nearest-neighbour distance),
so that a variable species forgives more. It also records the rank
where the two names part. Calls:

| Call | Meaning | Cases |
|---|---|---:|
| both | typical of both species: a real shared ITS2 (sisters, morphs) | 119 |
| fits one | typical of one species, unlike the other's records: the other label is the likely mislabel | 46 |
| one-sided | a name with no other records, nothing to judge by | 64 |
| neither | typical of neither | 4 |

220 of the 233 part at species. The 13 above species are the known
cross-genus suspects (C, §5.3).

**Rule for this group.** A record goes in only when all of these hold:
- the call is "fits one";
- the suspect label rests on one GenBank accession (no genome or long-read copy);
- both species have at least 2 other records to judge by;
- the pair does not also point the other way elsewhere;
- the depositing study is not itself about the two species.

Each row drops one source record. Where that record was a member of a
shared record, the record keeps the other name and becomes named.

| Accession | Deposited as | Fits | Split | Fits: records / nearest | Deposited as: records / nearest | Source |
|---|---|---|---|---:|---:|---|
| FJ534515 | *Meloidogyne incognita* | *Meloidogyne enterolobii* | species | 40 / 0 | 66 / 25 | unpublished, Malaysia |
| PV646252 | *Brugia malayi* | *Brugia pahangi* | species | 28 / 1 | 38 / 8 | unpublished, dog, Thailand; described as a novel zoonotic Brugia |
| PP953626 | *Contracaecum gibsoni* | *Contracaecum micropapillatum* | species | 7 / 0 | 2 / 5 | unpublished, pelican, Iran |
| EU180057 | *Hemicriconemoides californianus* | *Hemicriconemoides chitwoodi* | species | 6 / 0 | 4 / 80 | Chen et al., tea, Taiwan |
| KJ776346 | *Cloacina ernabella* | *Cloacina caenis* | species | 8 / 1 | 3 / 21 | Lott et al. 2014, T-RFLP study |
| KJ489433 | *Meloidogyne incognita* | *Meloidogyne enterolobii* | species | 40 / 0 | 66 / 25 | direct submission, India |
| KC608621 | *Steinernema kraussei* | *Steinernema silvaticum* | species | 2 / 1 | 41 / 21 | unpublished, Poland |
| PP953624 | *Contracaecum gibsoni* | *Contracaecum micropapillatum* | species | 7 / 1 | 2 / 5 | unpublished, pelican, Iran |
| KY088210 | *Heterorhabditis mexicana* | *Heterorhabditis sonorensis* | species | 8 / 1 | 6 / 4 | direct submission, Egypt |
| FN435301 | *Hemicycliophora conida* | *Hemicycliophora gracilis* | species | 4 / 1 | 3 / 60 | unpublished, a Mesocriconema xenoplax study |
| MN861048 | *Steinernema feltiae* | *Steinernema litorale* | species | 9 / 1 | 144 / 16 | direct submission, Turkey |
| OQ586460 | *Pratylenchus coffeae* | *Pratylenchus scribneri* | species | 30 / 0 | 145 / 115 | direct submission |
| MN341010 | *Heterorhabditis zealandica* | *Heterorhabditis indica* | species | 110 / 0 | 17 / 91 | direct submission |
| OP026096 | *Steinernema siamkayai* | *Steinernema abbasi* | species | 29 / 0 | 46 / 52 | unpublished, India |
| AF077090 | *Meloidogyne chitwoodi* | *Meloidogyne fallax* | species | 6 / 0 | 5 / 8 | Hugall et al., unpublished |
| AY045756 | *Heterodera hordecalis* | *Cactodera betulae* | genus | 2 / 1 | 18 / 62 | unpublished, Italy |
| JN808124 | *Steinernema intermedium* | *Steinernema affine* | species | 15 / 1 | 3 / 34 | unpublished, Portugal |
| MT193647 | *Trichostrongylus colubriformis* | *Trichostrongylus probolurus* | species | 4 / 2 | 47 / 6 | unpublished, Bactrian camel, China (a host of T. probolurus) |

Livestock and companion-animal rows: MT193647 (a camel "*T.
colubriformis*" that fits *T. probolurus*, a camel parasite) and PV646252
(a dog "*B. malayi*" that fits *B. pahangi*, the dog *Brugia*).

**Left out, for a person to judge:**
- ***Parascaris equorum* / *P. univalens*.** The triage calls *P.
  univalens* the mislabel, including genome-assembly copies. Most horse
  *Parascaris* sequenced today is probably *P. univalens*, and older
  "*P. equorum*" GenBank records may be too. A majority of labels is not
  evidence here.
- ***Ditylenchus dipsaci* / *D. weischeri*.** One record fits each way
  (MG386845, MG384729).
- ***Hysterothylacium fabri*.** Its own records are up to about 100 edits
  apart, so the label covers more than one lineage and "typical" means
  little (MW370925, MW370934, MW370888, MW370927).
- **MH030599** *Baylisascaris schroederi* (fits *B. transfuga*), deposited
  by a study of *B. transfuga* (Hoberg et al. 2018), so possibly a
  deliberate comparison. The two are close sisters (panda and bear).
- **AY702699** *Hovorkonema variegatum* (fits *Cyathostoma americana*,
  genus split) and **LN813018** *Trichuris globulosa* (fits *T. ovis*).
  Both come from studies of the specific status of exactly these taxa
  (Krone et al. 2007; Cutillas et al.). Their identity may be the
  studies' finding, not an error. Read the papers.
- ***Cylicocyclus insigne* / *C. ashworthi*.** Three long-read variants
  of one specimen (ERR11658755 v1, v4, v7) fit *C. ashworthi*. These are
  rDNA copies within one worm, which says the two species share ITS2
  types (§5.3). Not a mislabel.

**Effect.** Only 8 of the 43 exact-match errors in the leave-one-out
involve one of these pairs. The rest are "both": genuinely shared ITS2,
which the gap rule already stops at genus. This group tidies labels;
it does not change the picture much.

## D. Decisions

- ***Ascaris lumbricoides* / *A. suum*: kept apart** (decided
  2026-11-11). The records interleave at 0-1 edits, so reads stop at
  *Ascaris*. Reaffirmed after the GenBank update: 7 and 6 records
  (and one shared), 274-285 bases; between the species a median 2 edits
  (0-9), within them 1 (0-4) and 2 (0-8). An exact match to the two
  records furthest from the other species (*A. suum* PP179022, 6 edits;
  *A. lumbricoides* LN600404, 3) is still named.
- **NCBI look-alike pairs** (DESIGN §5.3): *Phocanema* / *Pseudoterranova*,
  *Labiosimplex* / *Labiostrongylus longispicularis*, *Gracilacus* /
  *Paratylenchus paralatescens*, *Japanema* / *Sheraphelenchus*. These
  are kept apart because NCBI keeps them apart. Overriding NCBI would
  remove their guaranteed "wrong species" calls.
- **The other 50-odd `sib=0` records**, mostly plant-parasitic and fish
  nematodes, are in `results/curation/outlier_records.txt`. Examples:
  - *Heterodera sonchophila* MW346595 holds 100 *H. glycines* records;
  - *Globodera pallida* DQ097514 holds 34 *G. rostochiensis* records.

  Whether to go through these depends on who uses the database.

## F. Records that do not sit with their own species (28 rows, applied)

Found by checking *Ashworthius sidemi* against NCBI. NC1 matches all
11 *A. sidemi* records that hold the 5.8S end, and NC2 the one that
holds its 28S site, so the primers are not the problem. But of the
three *A. sidemi* records here, the one held by 10 accessions from
independent studies is the species. OP787564 (European bison, Poland)
is 81 edits from it and 2 from *Haemonchus contortus*, and AJ577470
(unpublished) is 79 from it and 33+ from everything else.

`bin/species_coherence.py` asks the same of every named record in a
species with two or more records (`results/curation/species_coherence.tsv`):
`own`, the edits to the nearest other record of its species; `other`,
the edits to the nearest record of another species; and `to_core`, the
edits to the species' best-supported record (most accessions, >= 3).

- **displaced** (53): `other` <= 3 and `own` >= 10. 32 of them are one
  accession next to a species whose record sits with its own: the
  likely mislabels. 24 of the 32 are proposed for dropping. Held back:
  six *Hysterothylacium* records (MW370920, MW370925, MW370927,
  MW370746, MW370896, MW413285; a genus whose labels cover several
  lineages, as in E), *Ditylenchus dipsaci* MG384729 (as in E) and
  *Meloidogyne ethiopica* LN626967 (identical to *M. incognita*; the
  tropical root-knot species share near-identical ITS2). The other 21
  rest on two-record species or a neighbour that is itself odd: listed,
  not proposed.
- **stray** (171): `own` >= 20, >= 20 from the core, and nothing within
  10 edits. 104 sit in 15 species with three or more strays each
  (*Acanthocheilonema reconditum* 19, *Parabronema skrjabini* 16,
  *Angiostrongylus malaysiensis* 16): real diversity within the label,
  or cryptic species, not single mislabels. Left. Of the 67 isolated
  strays, four are proposed: AJ577470 (above), MK087177 (*H.
  contortus*, 119 bases, 39 edits from every other *H. contortus*
  record), and AJ242966 / AJ242967 (*Dirofilaria repens*: the 5S rRNA
  spacer, not this ITS2). The rest are listed. A stray holds no other
  species at genus; a read matching it would be named wrongly, but few
  reads will.
- **closer** (400): nearer another species than its own, by less. Sister
  species and partials; listed only.
- **layer** (23): long-read copies within one worm far from their
  species' other records (*Cylicocyclus nassatus* copies 17-88 edits
  off, *Strongylus vulgaris* 13-31). rDNA copies, not labels, but the
  far ones may be pseudogene copies; worth a look in the layer itself.

Also corrects C: dropping OP787564 does not remove *A. sidemi*.
The 28 rows are in `curation/overrides.tsv`.

## G. GenBank since the backbone (applied; 733 records in the I rebuild)

`bin/genbank_update.py` fetches nematode records the backbone (nemabiome
20260317) lacks: ITS2 records published since 17 March 2026 (1,773),
5.8S / 28S records published since then (1,919), and records of any date
annotated with 5.8S and 28S but never ITS2 (128). Of 3,820, 3,055 are
formally named species (not "sp.", cf., aff.), and 1,737 hold an ITS2
region (>= 80 bases at >= 70% to a record here). They are cut and graded
with the backbone (`bin/its2.py`, `bin/its2_partials.py`), which leaves
the backbone's own cuts as they were (15,103 of 15,108 identical).

`bin/its2.py` places its anchors near the backbone HMM's ITS2 boundary;
for a new record that boundary is where the record aligns to our records.
Where no anchor is found it falls back to that boundary, which for many
plant-parasitic and entomopathogenic records ran into ITS1 (*Heterodera*
median 539 bases against 212). So only the 601 records cut at an anchor
or at their own end are taken; the other 1,136 (mostly Tylenchida and
*Heterorhabditis*, and the 30 new *Dictyocaulus filaria*) wait for a
better boundary. Of the 601, 167 join existing records, 188 fold into
one, and 246 make 160 new records: 11,273 records, 1,869 species.
Ruminant gains: *Dictyocaulus viviparus* 19 -> 48 accessions (14 -> 43
records), *Cystocaulus ocreatus* 7 -> 37, *Oesophagostomum columbianum*
10 -> 21, *Haemonchus contortus* 476 -> 516, *Trichostrongylus
colubriformis* 155 -> 181, and a first complete *Ostertagia gruehneri*
(PZ782686).

Screened as the backbone was. 34 not to take in, now in
`curation/overrides.tsv`:
- 23 "*H. contortus*" from goats and sheep in Ghana (PZ4068xx, 2026)
  carrying the *H. placei* ITS2 of LC367244 and OR702556, which no *H.
  contortus* genome copy carries (DESIGN §5.3); taken in, they would make
  that *H. placei* record shared and hold its reads at genus;
- one *Steinernema kraussei* identical to *S. feltiae*;
- 10 strays (8 records), among them PZ369172 (*Ashworthius sidemi*,
  titled ITS1 + 5.8S).

## H. Records that still hold 5.8S (22, fixed in the cut)

22 records hold NC1 (<= 2 edits) more than 25 bases in, so 5.8S, and in
two cases ITS1, was left on: *Ostertagia antipini* KR259650 (elk; 656
bases, ITS1 + 5.8S + a few bases of ITS2), two *Heterodera*, eight
marine-mammal lungworms (*Stenurus*, *Parafilaroides*, 29 bases each),
seven *Caenorhabditis*, *Ancylostoma ceylanicum* AJ001593 (106 bases;
also displaced, F). `bin/its2.py` now cuts at an NC1 (<= 2 edits) more
than 25 bases inside a cut (DESIGN §3.2).

## I. The rebuild after v0.1.3 (applied)

The update rerun with `curation/genbank_named.tsv` (98 accessions found
missing by hand: strongylid species and variants, informal
*Oesophagostomum*, *Cloacina* and *Trichostrongylus* lineages; 1,834
records holding ITS2, 733 kept), cut with H's fix and a second one: a
record ending in a 28S stub too short for the anchor is cut at the stub
(DESIGN §3.2), which brings in complete *Oesophagostomum radiatum*
PZ347035/37. Screened again (`bin/species_coherence.py`; new records
within 2 edits of another species):

- **PZ569106** "*H. contortus*" (isolate sm001, Nigeria, 2026) carries
  the *H. placei* ITS2, 1 edit from the main *H. placei* ASVs, as the
  Ghana records in G. A single record, so the shared-record check missed
  it, and *Haemonchus* species sit 4 edits apart, so the coherence
  screen did too. It held PRJNA1524548's main *H. placei* ASVs at genus
  (2x150 reads to species 95.9% -> 71.2%). Dropped.
- **MK573551** "*Cloacina* sp. 2 AJ-2019" is identical to *Zoniolaimus
  occidentalis* MT084411 and OM987392, 28 edits from the nearest
  *Cloacina*: a *Zoniolaimus*. Kept, it pulled the *Z. occidentalis*
  record to family. Dropped.
- **MG711656** "*Trichostrongylus* sp. RR-2017" is 0 edits over its span
  from *Cooperia pectinata*, 23 from any *Trichostrongylus*. Dropped.
- **Informal names are not species.** `bin/build.py` read "sp." names
  as species (the backbone has none), so 19 informal records went in as
  named species of their own. "sp."/"spp." now count as no name, as cf.
  and hybrids do: the records are genus-level, `unnamed`, with the
  informal name in `lineage`.
- **Two pairs of informal *Oesophagostomum* isolates** have identical
  ITS2s (MW621346 / MW528477; PP235521 / PP235523, each pair from one
  study). Each pair shares one `lineage` name, joined with " / ".

Kept, as real neighbours: *Heligmosomoides bibullosus* 1 edit from *H.
vandegrifti*, *Metastrongylus asymmetricus* 1 from *M. salmi*, a complete
*T. retortaeformis* 4 from *T. axei* (3 through 2x150 windows), *T.
colubriformis* records 2 from *T. capricola* JF276023 (C).

**Partials** (`--drop-partials`, tested, not applied): dropping the 394
partial records (< 75% of the ITS2) of species with a full-length
record loses no species, and changes the gap rule's leave-one-out by
0.3 points of species recall or less (56.4% -> 56.7% whole) with
precision and over-classification unchanged. 277 of those partials are
>= 4 edits from their species' full-length records: variants a read can
match, which a test of full-length queries cannot count. Kept.

**Effect**, gap rule, v0.1.3 -> after I (`results/`):

| Test | v0.1.3 | After I |
|---|---:|---:|
| Leave-one-out, species present, right species (whole) | 55.3% | 56.4% |
| same, 2x150 | 58.1% | 59.3% |
| Leave-one-out, species present, wrong species (whole) | 0.7% | 0.4% |
| Leave-one-out, species missing, wrong species (whole) | 5.9% | 5.4% |
| same, 2x150 | 6.5% | 5.9% |
| Ruminant mock, species recall | 87.8% | 89.7% |
| Canine / equine mocks | 100 / 99.0% | same |
| Francis & Šlapeta 2x250 hybrid, reads to species | 93.8% | 93.9% |
| PRJNA1524548 2x250 hybrid, reads to species | 95.2% | 94.3% |

Missing-species errors (whole) fell from 187 to 173. assignTaxonomy and
IdTaxa on the same queries: species recall 86.0% to 86.3% and 78.8% to
80.3%, over-classification 43.5% to 42.9% and 23.2% to 20.6%.
PRJNA1524548 loses a point to two more honest calls: some *Oe. radiatum*
ASVs are 98.3-98.8% to the complete record that now holds the partials
they matched exactly, and *T. axei* through 2x150 windows is 3 edits from
the new *T. retortaeformis*.
