#!/usr/bin/env Rscript
# Merge per-cell salmon quant.sf into a gene-level matrix (genes x cells) with tximport.
# Usage: tximport.R tx2gene.tsv count|tpm <cell quant dirs...>
a <- commandArgs(TRUE)
unit <- a[2]
dirs <- sort(a[-(1:2)])
files <- setNames(file.path(dirs, "quant.sf"), basename(dirs))
tx <- tximport::tximport(files, type = "salmon", tx2gene = read.delim(a[1], header = FALSE), dropInfReps = TRUE)
m <- if (unit == "count") tx$counts else tx$abundance
write.table(m, paste0("gene_", unit, ".tsv"), sep = "\t", quote = FALSE, col.names = NA)
