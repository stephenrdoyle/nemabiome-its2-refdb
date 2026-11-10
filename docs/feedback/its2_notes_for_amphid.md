# Notes for amphid: using the nemabiome ITS2 database

From nemabiome-its2-refdb, scratch build of 2026-11-10. The
database follows vespa-18s-refdb's release format. These are the
places where ITS2 needs something V4 did not.

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
the species-level records, 410 are tied with another species, 900 are
one edit away and 711 two.

A leave-one-out test (DESIGN §6.1; 3,164 queries, 1,210 species) shows
identity alone is the wrong measure. Under the V4 limits only 18.5% of
reads whose species is in the database are named to species. What
separates a right call from a wrong one is the margin to the next
species. Recommended rule:

- **species**: >= 99% identity over the span AND the best reference of
  any other species is >= 4 edits further off (same span);
- **genus** >= 97%, **family** >= 95%;
- a reference with `species_status` = `unnamed` (hybrids, cf.,
  environmental) is not "another species"; a `shared` one is.

| Rule (whole reads) | Species missing: wrong species | Present: right species | Present: wrong species |
|---|---:|---:|---:|
| V4 limits | 2.4% | 18.5% | 1.4% |
| gap rule | 6.6% | 52.9% | 0.7% |
| dada2 assignTaxonomy, for scale | 42.4% | 85.1% | 4.6% |
| DECIPHER IdTaxa, for scale | 22.0% | 77.1% | 2.5% |

Limits set separately per genus or family did no better on held-out
species. Thresholds cannot fix the remaining missing-species errors:
they are reads identical to a record deposited under another name.

**Report how close a call came.** With each call, give the best species
and its identity, and the nearest other species with how many edits
further it is. Example: "*Trichostrongylus*; nearest *T. retortaeformis*
100%; next *T. axei* +2". A genus call then says whether it is "exactly
species X, but its sister is 2 edits away" or "98% to X, nothing else
close". Both stop at genus for different reasons, and users act on them
differently. `bin/classify_asvs.py` writes these as `nearest` and
`next`.

## 5a. If a k-mer classifier is used instead

assignTaxonomy and IdTaxa compare k-mers over the whole ASV, not over
the reference's span. Against this database, whose records stop before
the 28S stretch, every ASV carries ~47 bases no reference has. Cut the
stretch off the ASVs first (`bin/trim_28s.py`: at the 28S anchor
`AACCTGAGCTCAGTCGTGATTACCCGCCGAACTTAAGCATATCA`, else the last 47
bases). Uncut, assignTaxonomy named 45% of reads of one sheep set to
species; cut, 99.3%.

## 6. Ties are already resolved in the tables

A record whose span sits exactly inside another species' record is
"tied" with it. The resolution tables stop it at the shared rank. If
amphid assigns from the tables, nothing more is needed. If it assigns
from alignments, a read hitting both at the same score must be given
their common ancestor. That is the usual LCA rule, but it only works if
both hits are reported, so secondary alignments must not be capped
below the number of tied records.
