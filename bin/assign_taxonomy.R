#!/usr/bin/env Rscript
# dada2::assignTaxonomy(minBoot = 80) for an ASV table from bin/dada2_nemabiome.R,
# against one reference - to rerun a reference without rerunning DADA2.
#
#   Rscript bin/assign_taxonomy.R <asvs.tsv> <ref.fasta> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
lib <- file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))),
                 "work", "R-lib")
.libPaths(c(lib, .libPaths()))
suppressMessages(library(dada2))
set.seed(1)
a <- read.delim(args[1], stringsAsFactors = FALSE)
tx <- assignTaxonomy(a$sequence, args[2], minBoot = 80, outputBootstraps = TRUE, tryRC = FALSE,
                     multithread = as.integer(Sys.getenv("THREADS", "4")), verbose = FALSE)
write.table(data.frame(asv = a$asv, tx$tax, boot = tx$boot, check.names = FALSE), args[3],
            sep = "\t", quote = FALSE, row.names = FALSE)
