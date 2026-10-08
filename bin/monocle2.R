#!/usr/bin/env Rscript
# SIP downstream with Monocle 2: DDRTree trajectory and pseudotime.
# Usage: monocle2.R matrix.tsv genes.tsv count|tpm|fpkm
suppressMessages(library(monocle))
a <- commandArgs(TRUE)
set.seed(1)
m <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
g <- read.delim(a[2], header = FALSE, row.names = 1)
fd <- new("AnnotatedDataFrame", data.frame(
    gene_short_name = ifelse(rownames(m) %in% rownames(g), g[rownames(m), 1], rownames(m)), row.names = rownames(m)))
pd <- new("AnnotatedDataFrame", data.frame(cell = colnames(m), row.names = colnames(m)))

if (a[3] != "count") {   # Monocle 2's route for TPM/FPKM: convert to relative transcript counts (Census)
    m <- relative2abs(newCellDataSet(m, pd, fd, expressionFamily = tobit()), method = "num_genes")
}
cds <- newCellDataSet(round(m), pd, fd, expressionFamily = negbinomial.size())
cds <- estimateDispersions(estimateSizeFactors(cds))

disp <- dispersionTable(cds)
cds <- setOrderingFilter(cds, subset(disp, mean_expression >= 0.1 & dispersion_empirical >= dispersion_fit)$gene_id)
cds <- orderCells(reduceDimension(cds, max_components = 2, method = "DDRTree"))

pdf("trajectory.pdf")
print(plot_cell_trajectory(cds, color_by = "State"))
print(plot_cell_trajectory(cds, color_by = "Pseudotime") + guides(colour = guide_colourbar(display = "rectangles")))
dev.off()
xy <- t(reducedDimS(cds))   # DDRTree coordinates, for the summary report
write.table(data.frame(pData(cds)[, c("Pseudotime", "State")], Component1 = xy[, 1], Component2 = xy[, 2]),
            "pseudotime.tsv", sep = "\t", quote = FALSE, col.names = NA)
