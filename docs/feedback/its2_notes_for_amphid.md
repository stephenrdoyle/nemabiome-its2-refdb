# Notes for amphid: using the nemabiome ITS2 database

From nemabiome-its2-refdb v0.1.1 (2026-11-10). The
database follows vespa-18s-refdb's release format. These are the
places where ITS2 needs something V4 did not.

Tracked in amphid as
[stephenrdoyle/nf-amphid-pipeline#9](https://github.com/stephenrdoyle/nf-amphid-pipeline/issues/9),
which maps each point below to the amphid code it changes. The issue
also covers two points not here: `taxonomy.tsv` has NCBI's eight ranks,
not the PR2 levels amphid hard-codes, and amphid's identity charges a
read's unaligned ends - for ITS2, the 28S stretch - as differences.

## 1. Score reads over the reference's span

A reference runs from the 3' end of NC1 to the start of the ~47-base
28S stretch before NC2. Reads carry that stretch; references do not,
because it varies by 1-6 bases between nematodes and borrowing one
would invent differences. About two thirds of references are also
partial at one end or both (`grade` below).

So identity must be computed over the reference's span: the reference
aligned end to end, the read's bases beyond it free (semi-global /
"glocal" alignment, reference global, read local). Scoring over the
read's full length would charge every read ~47 bases of 28S against
every reference, and charge a partial reference for bases it simply
lacks.

`manifest.json` states this as `score_over`. Margins and resolution
tables in the release are computed the same way.

## 2. Primers

NC1 `ACGTCTGGTTCAGGGTTGTT`, NC2 `TTAGTTTCTTTTCCTCCGCT`, each ordered
with 0-3 N spacers behind the Illumina adapter (nemabiome.ca's
protocol). Primer trimming must allow a 0-3 base offset.

## 3. Read windows

2x150, primers trimmed: forward read = first 130 bases of the ITS2;
reverse read = last 130 bases of the amplicon, of which ~47 are 28S, so
the last 83 of the ITS2. 2x250: 230 + 183. With a median ITS2 of 240
bases, most 2x250 pairs and many 2x150 pairs overlap and merge.
`manifest.json` → `illumina_geometries` gives the windows.

## 4. Grade

`taxonomy.tsv` and the resolution tables carry each record's `grade`:
`complete`, `near_complete` (>= 98% of the ITS2), `partial_90`,
`partial_75`, `partial` (< 75%), `unmeasured`. A read matching a
`partial` reference at 100% has been checked on fewer bases. Report the
grade with the call, and consider requiring a minimum aligned length
(say 150 bases) for a species call.

## 5. Identity limits: use the gap rule

The V4 limits (species 99.7% Illumina / 99.5% nanopore, genus 99.5% /
99%, family 95%) were set for a ~400-base insert. Over a 240-base ITS2,
one edit is 0.42%, so 99.7% means an exact match. Margins are thin: of
the species-level records, 410 are tied with another species, 897 are
one edit away and 712 two.

A leave-one-out test (DESIGN §6.1; 3,182 queries, 1,220 species, on
v0.1.2) shows identity alone is the wrong
measure. Under the V4 limits only 18.2% of reads whose species is in
the database are named to species. What separates a right call from a
wrong one is the margin to the next species. Recommended rule:

- **species**: >= 99% identity over the span AND the best reference of
  any other species is >= 4 edits further off (same span); or an exact
  match (0 edits) AND that reference is >= 2 edits off;
- **genus** >= 97%, **family** >= 95%;
- a reference with `species_status` = `unnamed` (hybrids, cf.,
  environmental, informal "sp." names) is not "another species"; a
  `shared` one is.

Recall = right calls / reads; precision = right calls / calls made at
that rank; over-classification = a species named for a read whose
species is missing from the database.

| Rule (whole reads) | Species recall | Species precision | Genus recall | Genus precision | Over-classification |
|---|---:|---:|---:|---:|---:|
| V4 limits | 18.5% | 95.2% | 50.5% | 99.3% | 1.8% |
| gap rule | 56.4% | 99.2% | 89.4% | 99.4% | 5.4% |
| dada2 assignTaxonomy, for scale | 86.3% | 95.4% | 97.7% | 99.3% | 42.9% |
| DECIPHER IdTaxa, for scale | 80.3% | 97.5% | 94.5% | 99.6% | 20.6% |

(`main` after v0.1.3; v0.1.3 gave the gap rule 55.3%, 98.8% and 5.9%.)

Limits set separately per genus or family did no better on held-out
species. Thresholds cannot fix the remaining missing-species errors:
four in five are sister species under 1% apart with one of them
missing, the rest reads identical to a record deposited under another
name.

**Report how close a call came.** With each call, give the best species
and its identity, and the nearest other species with how many edits
further it is. Example: "*Trichostrongylus*; nearest *T. colubriformis* 100%;
next *T. capricola* +1". A genus call then says whether it is "exactly
species X, but its sister is 1 edit away" or "98% to X, nothing else
close". Both stop at genus for different reasons, and users act on them
differently. `bin/classify_asvs.py` writes these as `nearest` and
`next`.

**Report the lineage.** `taxonomy.tsv` has a `lineage` column: a
published line within a species, kept under the species' name: the
*Teladorsagia circumcincta* goat line (2 records), which Francis &
Šlapeta (2023) found in goats and not in sheep on shared pasture. On
`main` after v0.1.3 it also carries the informal name of a genus-level
record (19: *Oesophagostomum* sp. 3335_MDR_CAM_OH15, *Cloacina* sp. 1
AJ-2019, *Trichostrongylus* sp. SMNPOL33R), whose species is blank. When every tied best reference carries the same lineage,
report it with the call, at any rank: a goat-line read at 98% stops at
genus but is still "nearest *T. circumcincta* (goat line) 98.0%". Blank
means no line is recorded, not the other line. `bin/classify_asvs.py`
writes it as `lineage`.

**Score read windows strictly.** Charge every reference base from the
reference's own end to where a window stops; leave read bases free only
before the reference's start (primer offset, up to 30). Cut the second
window at the 28S anchor, not at a fixed 47 bases. Giving a reference's
first bases slack let one partial *H. bedfordi* record tie *H. placei*
reads in the reference scorer (DESIGN §6.2); amphid already charges
them.

**Orient reads.** A sequence whose best hit is on the reverse strand
should be turned before scoring; IdTaxa does this with
`strand = "both"`.

## 5a. If a k-mer classifier is used instead

assignTaxonomy and IdTaxa compare k-mers over the whole ASV, not over
the reference's span. Against this database, whose records stop before
the 28S stretch, every ASV carries ~47 bases no reference has. Cut the
stretch off the ASVs first (`bin/trim_28s.py`, or `R/trim_28s.R` as
shipped with the release: at the 28S anchor
`AACCTGAGCTCAGTCGTGATTACCCGCCGAACTTAAGCATATCA`, else the last 47
bases). Uncut, assignTaxonomy named 45% of reads of one sheep set to
species; cut, 99.2%.

## 6. Ties are already resolved in the tables

A record whose span sits exactly inside another species' record is
"tied" with it. The resolution tables stop it at the shared rank. If
amphid assigns from the tables, nothing more is needed. If it assigns
from alignments, a read hitting both at the same score must be given
their common ancestor. That is the usual LCA rule, but it only works if
both hits are reported, so secondary alignments must not be capped
below the number of tied records.
