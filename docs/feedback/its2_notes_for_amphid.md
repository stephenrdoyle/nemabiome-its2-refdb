# Notes for amphid: using the nemabiome ITS2 database

From nemabiome-its2-refdb, first scratch build (2026-11-09). The
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

## 5. Identity limits

The V4 limits (species 99.7% Illumina / 99.5% nanopore, genus 99.5% /
99%, family 95%) were set for a ~400-base insert. Over a 240-base ITS2,
one edit is 0.42%, so 99.7% means an exact match and 99.5% means at
most one edit. Margins are thin: of 12,857 species-level records, 732
are tied with another species (margin 0), 1,040 are one edit away and
859 two.

Two suggestions:

- express the species limit in edits over the aligned span (0 for
  Illumina, <= 1 for nanopore after correction), or keep percent
  but check it against the margin table;
- a call whose best and second-best references (different species)
  are within one edit should stop at genus, whatever the identity.
  `species_margins.tsv` gives each record's nearest other species.

A leave-one-out test (as for V4) is planned and will set these
properly.

## 6. Ties are already resolved in the tables

A record whose span sits exactly inside another species' record is
"tied" with it. The resolution tables stop it at the shared rank. If
amphid assigns from the tables, nothing more is needed. If it assigns
from alignments, a read hitting both at the same score must be given
their common ancestor. That is the usual LCA rule, but it only works if
both hits are reported, so secondary alignments must not be capped
below the number of tied records.
