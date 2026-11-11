#!/usr/bin/env Rscript
# The standard nemabiome DADA2 workflow (nemabiome.ca, after Callahan et
# al. 2016), run on primer-trimmed read pairs, for comparison with this
# database's span-scored calls.
#
#   Rscript bin/dada2_nemabiome.R <trimmed fastq dir> <out prefix> <modes> [<ref.fasta> ...]
#
# <modes>, one or more, comma-separated; denoising is shared:
# merge:  mergePairs as usual (2x250). Pairs that do not overlap by 12
#         bases are lost: an ITS2 over ~400 bases never merges at 2x250
#         (Metastrongylus, Trichuris, filarial worms).
# concat: mergePairs(justConcatenate = TRUE) - forward + 10 N + reverse -
#         for every pair (2x150 over a ~290-330 base amplicon, or 2x250
#         without merging). The number that would have merged is reported.
# hybrid: merged where the pair overlaps, concatenated where it does not.
# With more than one mode, each writes <out prefix>.<mode>.*.
#
# filterAndTrim(maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50),
# learnErrors, dada, merge, removeBimeraDenovo(consensus), then
# assignTaxonomy(minBoot = 80) against each reference given.
#
# Writes <prefix>.asvs.tsv (asv id, sequence), <prefix>.counts.tsv
# (sample, asv, reads), <prefix>.track.tsv (reads kept at each step) and
# <prefix>.tax.<ref basename>.tsv (asv, ranks, bootstraps).

args <- commandArgs(trailingOnly = TRUE)
lib <- file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))))),
                 "work", "R-lib")
.libPaths(c(lib, .libPaths()))
suppressMessages(library(dada2))

dir <- args[1]; prefix <- args[2]; mode <- args[3]; refs <- args[-(1:3)]
threads <- as.integer(Sys.getenv("THREADS", "6"))
set.seed(1)

f <- sort(list.files(dir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
r <- sub("_1\\.fastq\\.gz$", "_2.fastq.gz", f)
samples <- sub("_1\\.fastq\\.gz$", "", basename(f))
ff <- file.path(dir, "filtered", basename(f)); fr <- file.path(dir, "filtered", basename(r))
out <- filterAndTrim(f, ff, r, fr, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50,
                     rm.phix = TRUE, compress = TRUE, multithread = threads)
keep <- out[, 2] > 0
ff <- ff[keep]; fr <- fr[keep]; samples <- samples[keep]

errF <- learnErrors(ff, multithread = threads, verbose = 0)
errR <- learnErrors(fr, multithread = threads, verbose = 0)
dF <- dada(ff, err = errF, multithread = threads, verbose = 0)
dR <- dada(fr, err = errR, multithread = threads, verbose = 0)
if (length(ff) == 1) { dF <- list(dF); dR <- list(dR) }

one <- function(x) if (length(ff) == 1) list(x) else x
merged <- one(mergePairs(dF, ff, dR, fr, returnRejects = TRUE, verbose = FALSE))
concat <- one(mergePairs(dF, ff, dR, fr, justConcatenate = TRUE, verbose = FALSE))
n_merged <- sapply(merged, function(m) sum(m$abundance[m$accept]))
denoised <- sapply(dF, function(d) sum(getUniques(d)))
modes <- strsplit(mode, ",")[[1]]

write_mode <- function(m) {
  pairs <- lapply(seq_along(merged), function(i) {
    acc <- merged[[i]][merged[[i]]$accept, ]
    if (m == "merge") return(acc)
    if (m == "concat") return(concat[[i]])
    key <- paste(concat[[i]]$forward, concat[[i]]$reverse)
    rbind(acc[, c("sequence", "abundance")],
          concat[[i]][!key %in% paste(acc$forward, acc$reverse), c("sequence", "abundance")])
  })
  names(pairs) <- samples
  st <- makeSequenceTable(pairs)
  nochim <- removeBimeraDenovo(st, method = "consensus", multithread = threads, verbose = FALSE)
  pre <- if (length(modes) > 1) paste0(prefix, ".", m) else prefix
  track <- data.frame(sample = samples, input = out[keep, 1], filtered = out[keep, 2],
                      denoised = denoised, merged_overlap = n_merged,
                      used = rowSums(st), nonchim = rowSums(nochim))
  write.table(track, paste0(pre, ".track.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  asv <- colnames(nochim)
  ids <- sprintf("ASV%05d", seq_along(asv))
  write.table(data.frame(asv = ids, sequence = asv), paste0(pre, ".asvs.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cnt <- which(nochim > 0, arr.ind = TRUE)
  write.table(data.frame(sample = rownames(nochim)[cnt[, 1]], asv = ids[cnt[, 2]], reads = nochim[cnt]),
              paste0(pre, ".counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  for (ref in refs) {
    tx <- assignTaxonomy(asv, ref, minBoot = 80, outputBootstraps = TRUE,
                         tryRC = FALSE, multithread = threads, verbose = FALSE)
    tab <- data.frame(asv = ids, tx$tax, boot = tx$boot, check.names = FALSE)
    write.table(tab, paste0(pre, ".tax.", sub("\\.fasta$", "", basename(ref)), ".tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  cat(sprintf("%s: %d samples, %d ASVs, %s reads merged by overlap, %s reads used\n", pre,
              length(samples), length(asv), format(sum(n_merged), big.mark = ","),
              format(sum(nochim), big.mark = ",")))
}
for (m in modes) write_mode(m)
