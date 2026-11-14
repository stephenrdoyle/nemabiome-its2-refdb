# trim_28s(): cut the 28S stretch off NC1/NC2 ASVs before classifying
# them with assignTaxonomy or IdTaxa against nematode-its2-refdb.
#
# Why: a primer-trimmed NC1/NC2 ASV ends in ~47 bases of 28S before NC2.
# This database's records stop where that stretch starts (it varies by
# 1-6 bases between nematodes, and most GenBank records lack it). k-mer
# classifiers compare the whole ASV, so uncut ASVs carry ~47 bases no
# reference has and lose confidence: on one sheep data set, assignTaxonomy
# named 45% of reads to species uncut and 99.3% cut.
#
# Usage, after removeBimeraDenovo:
#
#   source("trim_28s.R")
#   seqs <- getSequences(seqtab.nochim)
#   cut  <- trim_28s(seqs)
#   ok   <- nchar(gsub("N", "", cut)) >= 100      # shorter is not ITS2
#   # assignTaxonomy
#   taxa <- assignTaxonomy(cut[ok], "nemits2_dada2.fasta", minBoot = 80)
#   rownames(taxa) <- seqs[ok]                      # back to the ASVs
#   # or IdTaxa
#   load("nemits2_idtaxa.RData")                    # trainingSet
#   ids <- IdTaxa(DNAStringSet(cut[ok]), trainingSet, strand = "both",
#                 threshold = 60)
#
# The cut is at the 28S anchor AACCTGAGCTCAGTCGTGATTACCCGCCGAACTTAAGCAT,
# which stops at the end of the conserved motif so that reads cut there
# still fit it in place (<= 12 edits, best fit within the last 90
# bases), taken only if it removes 35-60 bases; otherwise the last 47
# bases go (anchor_only = TRUE leaves such an ASV whole). A concatenated
# ASV (forward + 10 N + reverse, from mergePairs(justConcatenate = TRUE))
# is cut in its reverse part.
# The same rule as bin/trim_28s.py.

S28_ANCHOR <- "AACCTGAGCTCAGTCGTGATTACCCGCCGAACTTAAGCAT"

# Start (0-based) of the best fit of pattern anywhere in text, the
# pattern aligned end to end; the earliest end on a tie. NA if no fit
# within max_edits. An N in the text matches any base.
.anchor_start <- function(text, pattern, max_edits) {
  p <- strsplit(pattern, "")[[1]]
  t <- strsplit(text, "")[[1]]
  m <- length(p)
  prev <- 0:m
  best <- max_edits + 1L
  best_end <- -1L
  for (j in seq_along(t)) {
    cur <- integer(m + 1)
    tj <- t[j]
    for (i in seq_len(m)) {
      cost <- if (tj == p[i] || tj == "N") 0L else 1L
      cur[i + 1] <- min(prev[i] + cost, prev[i + 1] + 1L, cur[i] + 1L)
    }
    if (cur[m + 1] < best) {
      best <- cur[m + 1]
      best_end <- j
    }
    prev <- cur
  }
  if (best_end < 0 || best > max_edits) return(NA_integer_)
  max(0L, best_end - m)
}

trim_28s <- function(seqs, anchor_only = FALSE, look = 90L, max_edits = 12L,
                     cut_range = c(35L, 60L), stretch = 47L) {
  sep <- strrep("N", 10)
  out <- vapply(toupper(seqs), function(s) {
    at <- regexpr(sep, s, fixed = TRUE)
    if (at > 0) {
      head <- substr(s, 1, at + nchar(sep) - 1)
      part <- substr(s, at + nchar(sep), nchar(s))
    } else {
      head <- ""
      part <- s
    }
    n <- nchar(part)
    cut <- NA_integer_
    if (n >= look) {
      st <- .anchor_start(substr(part, n - look + 1, n), S28_ANCHOR, max_edits)
      if (!is.na(st)) cut <- n - look + st
    }
    found <- !is.na(cut) && (n - cut) >= cut_range[1] && (n - cut) <= cut_range[2]
    if (!found) {
      if (anchor_only) return(s)
      cut <- n - stretch
    }
    paste0(head, substr(part, 1, cut))
  }, character(1), USE.NAMES = FALSE)
  names(out) <- names(seqs)
  out
}
