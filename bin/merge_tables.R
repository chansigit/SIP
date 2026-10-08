#!/usr/bin/env Rscript
# Merge per-cell featureCounts / HTSeq / StringTie tables into a gene-level matrix (genes x cells).
# Usage: merge_tables.R featurecounts|htseq|stringtie count|tpm|fpkm <tables...>
a <- commandArgs(TRUE)
tool <- a[1]; unit <- a[2]; files <- sort(a[-(1:2)])

read1 <- function(f) switch(tool,
    featurecounts = { d <- read.delim(f, comment.char = "#"); setNames(d[[7]], d$Geneid) },
    htseq = { d <- read.delim(f, header = FALSE); d <- d[!startsWith(d$V1, "__"), ]; setNames(d$V2, d$V1) },
    stringtie = { d <- read.delim(f, check.names = FALSE); tapply(d[[toupper(unit)]], d[["Gene ID"]], sum) })

cols <- lapply(files, read1)
genes <- sort(unique(unlist(lapply(cols, names))))
m <- sapply(cols, function(v) { x <- v[genes]; x[is.na(x)] <- 0; x })
dimnames(m) <- list(genes, sub("\\.(featureCounts\\.txt|htseq\\.txt|gene_abund\\.tab)$", "", basename(files)))
write.table(m, paste0("gene_", unit, ".tsv"), sep = "\t", quote = FALSE, col.names = NA)

if (tool == "stringtie") {   # merged-assembly gene ids (MSTRG.*) -> reference gene names
    d <- do.call(rbind, lapply(files, function(f) read.delim(f, check.names = FALSE)[, 1:2]))
    d <- d[d[[2]] != "-" & !duplicated(d[[1]]), ]
    write.table(d, "genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
