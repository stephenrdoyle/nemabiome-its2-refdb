#!/bin/bash
# The other arms for one simulated mock set (docs/DESIGN.md §6.3), after
# bin/dada2_nemabiome.R has written work/mocks/<name>/<name>.* (with
# assignTaxonomy against both references), then bin/score_mocks.py.
#
#   bin/mock_arms.sh <name> <source fasta in work/mocks/Multifasta_for_ART_Illumina>
#
# ASVs are cut at the 28S anchor only where it is found (--anchor-only):
# the mocks were simulated from database records, not all reaching 28S.
set -e
n=$1; P=work/mocks/$n/$n
python3 -I bin/trim_28s.py $P.asvs.tsv $P.asvs.trim28s.tsv --anchor-only
THREADS=4 Rscript bin/assign_taxonomy.R $P.asvs.trim28s.tsv work/mocks/nemits2_dev.fasta $P.tax.nemits2_dev.trim28s.tsv 2>/dev/null
THREADS=4 Rscript bin/idtaxa.R $P.asvs.tsv work/test/refs/nemabiome_20260317.fasta $P.idtaxa.nemabiome_20260317.tsv 2>/dev/null
THREADS=4 Rscript bin/idtaxa.R $P.asvs.trim28s.tsv work/mocks/nemits2_dev.fasta $P.idtaxa.nemits2_dev.trim28s.tsv 2>/dev/null
python3 -I bin/classify_asvs.py $P.asvs.tsv work/build $P.ours.nemits2_dev.tsv --workers 4 2>/dev/null
python3 -I bin/score_mocks.py work/mocks/$n work/mocks/Multifasta_for_ART_Illumina/$2 curation/overrides.tsv sources/nemabiome_its2_20260317/Nematode_ITS2_20260317_seqs_nr.fasta
