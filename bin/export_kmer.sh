#!/bin/bash
# The k-mer classifier files of a build, for the DADA2 workflow
# (README: For DADA2 and IdTaxa users):
#   <out>/<name>_dada2.fasta    assignTaxonomy reference (bin/export_dada2.py)
#   <out>/<name>_idtaxa.RData   IdTaxa trainingSet (bin/train_idtaxa.R)
#   <out>/trim_28s.R            cut the 28S stretch off ASVs first (required)
#
#   bin/export_kmer.sh <build dir> <out dir> <name>
set -euo pipefail
build=$1; out=$2; name=$3
here=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$out"
python3 -I "$here/bin/export_dada2.py" "$build" "$out/${name}_dada2.fasta"
Rscript "$here/bin/train_idtaxa.R" "$out/${name}_dada2.fasta" "$out/${name}_idtaxa.RData"
cp "$here/R/trim_28s.R" "$out/trim_28s.R"
