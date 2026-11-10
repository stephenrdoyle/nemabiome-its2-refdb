#!/usr/bin/env Rscript
# DECIPHER IdTaxa as the nemabiome.ca DADA2 workflow runs it: LearnTaxa on
# the reference, then IdTaxa(strand = "both", threshold = 60,
# bootstraps = 100). Takes a dada2-format reference (headers of
# ;-separated ranks, as bin/export_dada2.py and nemabiome's dada2 export
# write) and adds the "Root;" LearnTaxa needs.
#
#   Rscript bin/idtaxa.R <asvs.tsv> <ref.fasta> <out.tsv>
#
# Output as bin/assign_taxonomy.R's: asv, one column per rank (NA below
# the confidence threshold), then the confidences.
args <- commandArgs(trailingOnly = TRUE)
lib <- file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))),
                 "work", "R-lib")
.libPaths(c(lib, .libPaths()))
suppressMessages(library(DECIPHER))
set.seed(1)
ranks <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")

train <- readDNAStringSet(args[2])
tax <- paste0("Root;", sub(";$", "", names(train)))
train <- train[!duplicated(paste(tax, as.character(train)))]
tax <- paste0("Root;", sub(";$", "", names(train)))
ts <- LearnTaxa(train, tax, verbose = FALSE)

a <- read.delim(args[1], stringsAsFactors = FALSE)
dna <- DNAStringSet(gsub("[^ACGT]", "N", toupper(a$sequence)))
ids <- IdTaxa(dna, ts, strand = "both", threshold = 60, bootstraps = 100,
              processors = as.integer(Sys.getenv("THREADS", "4")), verbose = FALSE, type = "extended")
rows <- lapply(ids, function(x) {
  taxa <- x$taxon[-1]; conf <- x$confidence[-1]
  keep <- !grepl("^unclassified", taxa)
  taxa <- taxa[keep]; conf <- conf[keep]
  length(taxa) <- length(ranks); length(conf) <- length(ranks)
  c(taxa, round(conf, 1))
})
m <- do.call(rbind, rows)
out <- data.frame(asv = a$asv, m, check.names = FALSE, stringsAsFactors = FALSE)
names(out) <- c("asv", ranks, paste0("boot.", ranks))
write.table(out, args[3], sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
