#!/bin/bash
# The real-data comparison (docs/DESIGN.md §6.2), raw reads to arms, for
# one project:
#
#   bin/real_data.sh <project> <build dir> <out dir>
#
# Raw pairs in work/testdata/<project>/ (2x251). Primer-trimmed at 2x250
# and, truncated first, at 2x150 (bin/trim_primers.py). DADA2
# (bin/dada2_nemabiome.R): 2x150 concatenated; 2x250 three ways, merged
# (the standard workflow), concatenated, and hybrid (merged where the
# pair overlaps). Then on each ASV table: assignTaxonomy against the
# backbone and this build (in the DADA2 run), IdTaxa against both, both
# again on 28S-cut ASVs against this build, bin/classify_asvs.py against
# both, and bin/compare_arms.py. Steps whose output exists are skipped:
# delete an output to redo it.
set -e
P=$1; BUILD=$2; OUT=$3
T=${THREADS:-8}
REFS=$OUT/refs
mkdir -p $OUT $REFS
[ -s $REFS/nemits2_dev.fasta ] || python3 -I bin/export_dada2.py $BUILD $REFS/nemits2_dev.fasta
[ -s $REFS/nemabiome_20260317.fasta ] || cp work/test/refs/nemabiome_20260317.fasta $REFS/

[ -s $OUT/$P/250/trim_stats.tsv ] || python3 -I bin/trim_primers.py work/testdata/$P $OUT/$P/250 --workers $T
[ -s $OUT/$P/150/trim_stats.tsv ] || python3 -I bin/trim_primers.py work/testdata/$P $OUT/$P/150 --length 150 --workers $T

arms() {
  local X=$1
  [ -s $X.asvs.trim28s.tsv ] || python3 -I bin/trim_28s.py $X.asvs.tsv $X.asvs.trim28s.tsv
  [ -s $X.tax.nemits2_dev.trim28s.tsv ] || THREADS=$T Rscript bin/assign_taxonomy.R $X.asvs.trim28s.tsv $REFS/nemits2_dev.fasta $X.tax.nemits2_dev.trim28s.tsv 2>/dev/null
  [ -s $X.idtaxa.nemabiome_20260317.tsv ] || THREADS=$T Rscript bin/idtaxa.R $X.asvs.tsv $REFS/nemabiome_20260317.fasta $X.idtaxa.nemabiome_20260317.tsv 2>/dev/null
  [ -s $X.idtaxa.nemits2_dev.tsv ] || THREADS=$T Rscript bin/idtaxa.R $X.asvs.tsv $REFS/nemits2_dev.fasta $X.idtaxa.nemits2_dev.tsv 2>/dev/null
  [ -s $X.idtaxa.nemits2_dev.trim28s.tsv ] || THREADS=$T Rscript bin/idtaxa.R $X.asvs.trim28s.tsv $REFS/nemits2_dev.fasta $X.idtaxa.nemits2_dev.trim28s.tsv 2>/dev/null
  [ -s $X.ours.nemits2_dev.tsv ] || python3 -I bin/classify_asvs.py $X.asvs.tsv $BUILD $X.ours.nemits2_dev.tsv --workers $T 2>/dev/null
  [ -s $X.ours.nemabiome_20260317.tsv ] || python3 -I bin/classify_asvs.py $X.asvs.tsv $REFS/nemabiome_20260317.fasta $X.ours.nemabiome_20260317.tsv --workers $T 2>/dev/null
  python3 -I bin/compare_arms.py $X > $X.compare.txt
  echo "done $X"
}

[ -s $OUT/$P.150.counts.tsv ] || THREADS=$T Rscript bin/dada2_nemabiome.R $OUT/$P/150 $OUT/$P.150 concat \
  $REFS/nemabiome_20260317.fasta $REFS/nemits2_dev.fasta
arms $OUT/$P.150
[ -s $OUT/$P.250.hybrid.counts.tsv ] || THREADS=$T Rscript bin/dada2_nemabiome.R $OUT/$P/250 $OUT/$P.250 merge,concat,hybrid \
  $REFS/nemabiome_20260317.fasta $REFS/nemits2_dev.fasta
for m in merge concat hybrid; do arms $OUT/$P.250.$m; done
