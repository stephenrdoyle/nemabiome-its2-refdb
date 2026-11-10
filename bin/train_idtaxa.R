#!/usr/bin/env Rscript
# Train a DECIPHER IdTaxa classifier on a dada2-format reference
# (bin/export_dada2.py), for release beside it: users load() the file
# and get `trainingSet`, as for nemabiome.ca's own IdTaxa set.
#
#   Rscript bin/train_idtaxa.R <ref dada2.fasta> <out .RData>
#
# Identical sequence + taxonomy pairs are trained once; "Root;" is added
# as LearnTaxa needs. The same training as bin/idtaxa.R.
args <- commandArgs(trailingOnly = TRUE)
lib <- file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))),
                 "work", "R-lib")
.libPaths(c(lib, .libPaths()))
suppressMessages(library(DECIPHER))
set.seed(1)
train <- readDNAStringSet(args[1])
tax <- paste0("Root;", sub(";$", "", names(train)))
keep <- !duplicated(paste(tax, as.character(train)))
trainingSet <- LearnTaxa(train[keep], tax[keep], verbose = FALSE)
save(trainingSet, file = args[2], compress = "xz")
cat(sprintf("%s: %d sequences, DECIPHER %s\n", args[2], sum(keep), packageVersion("DECIPHER")))
