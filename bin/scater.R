#!/usr/bin/env Rscript
# SIP downstream with scater: per-cell QC metrics and overview plots.
# Usage: scater.R matrix.tsv genes.tsv
suppressMessages(library(scater))
a <- commandArgs(TRUE)
set.seed(1)
m <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
g <- read.delim(a[2], header = FALSE, row.names = 1)
rownames(m) <- make.unique(ifelse(rownames(m) %in% rownames(g), g[rownames(m), 1], rownames(m)))

sce <- SingleCellExperiment(list(counts = m))
sce <- addPerCellQCMetrics(sce, subsets = list(Mito = grep("^MT-", rownames(sce), ignore.case = TRUE)))
sce <- runPCA(logNormCounts(sce), ncomponents = min(10, ncol(sce) - 1))

pdf("scater_qc.pdf")
print(plotHighestExprs(sce))
print(plotColData(sce, x = "sum", y = "detected"))
print(plotPCA(sce, colour_by = "detected") + guides(colour = guide_colourbar(display = "rectangles"), fill = guide_colourbar(display = "rectangles")))
dev.off()
write.table(as.data.frame(colData(sce)), "cell_qc_metrics.tsv", sep = "\t", quote = FALSE, col.names = NA)
