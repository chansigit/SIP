#!/usr/bin/env Rscript
# SIP summary report: one vector PDF with the run settings, QC tables and the downstream figures.
# Usage: report.R 'key=value' ...   (run settings; result files are picked up by name from the working directory)
suppressMessages({ library(Seurat); library(ggplot2); library(gridExtra); library(grid) })
a <- commandArgs(TRUE)
rd <- function(f) read.delim(f, check.names = FALSE)
bar <- guide_colourbar(display = "rectangles")
title <- function(x) textGrob(x, gp = gpar(fontsize = 16, fontface = "bold"))
wrap <- function(x) vapply(x, function(s) paste(strwrap(gsub("[-_]", " ", s), 14), collapse = "\n"), "")
sig <- function(df) { df[] <- lapply(df, function(x) if (is.numeric(x)) signif(x, 4) else x); df }

table_pages <- function(head, df, rows = 28, cols = 9) {   # long or wide tables continue on the next page
    for (j in seq(2, max(ncol(df), 2), by = cols - 1)) for (i in seq(1, max(nrow(df), 1), by = rows)) {
        part <- df[i:min(i + rows - 1, nrow(df)), unique(c(1, j:min(j + cols - 2, ncol(df)))), drop = FALSE]
        names(part) <- wrap(names(part))
        grid.newpage()
        grid.draw(arrangeGrob(tableGrob(part, rows = NULL, theme = ttheme_minimal(base_size = 9)), top = title(head)))
    }
}

pdf("sip_report.pdf", width = 11, height = 8.5)

# run settings
mfile <- list.files(pattern = "^gene_.*\\.tsv$")[1]
m <- read.delim(mfile, row.names = 1, check.names = FALSE)
settings <- data.frame(item = c("date", sub("=.*", "", a), "expression matrix"),
                       value = c(format(Sys.time(), "%Y-%m-%d %H:%M"), sub("^[^=]*=", "", a),
                                 sprintf("%s: %d genes x %d cells", mfile, nrow(m), ncol(m))))
table_pages("SIP analysis report", settings)

# read and alignment QC
gs <- "multiqc_data/multiqc_general_stats.txt"
if (file.exists(gs)) table_pages("Read and alignment QC (MultiQC general statistics)", sig(rd(gs)))

if (file.exists("infer_experiment.txt")) {
    grid.newpage()
    grid.draw(arrangeGrob(textGrob(paste(readLines("infer_experiment.txt"), collapse = "\n"),
                                   gp = gpar(fontfamily = "mono", fontsize = 11)),
                          top = title("Strandness (RSeQC infer_experiment.py)")))
}

# Seurat
if (file.exists("seurat_plots.rds")) {
    cl <- rd("clusters.tsv")
    table_pages("Seurat: cells per cluster", as.data.frame(table(cluster = cl$cluster), responseName = "cells"))
    if (file.exists("markers.tsv")) {
        mk <- rd("markers.tsv")
        top <- do.call(rbind, lapply(split(mk, mk$cluster), function(d) head(d[order(-d$avg_log2FC), ], 10)))
        table_pages("Seurat: top 10 marker genes per cluster",
                    sig(top[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]))
    }
    figs <- readRDS("seurat_plots.rds")
    for (n in names(figs)) print(figs[[n]] + patchwork::plot_annotation(title = n))
}

# Monocle 2
if (file.exists("pseudotime.tsv")) {
    d <- rd("pseudotime.tsv")
    print(ggplot(d, aes(Component1, Component2, colour = Pseudotime)) + geom_point(size = 3) + guides(colour = bar) +
          theme_classic() + labs(title = "Monocle 2: DDRTree trajectory"))
    table_pages("Monocle 2: cells per state", as.data.frame(table(state = d$State), responseName = "cells"))
}

# scater
if (file.exists("cell_qc_metrics.tsv")) {
    d <- rd("cell_qc_metrics.tsv")
    print(ggplot(d, aes(sum, detected, colour = subsets_Mito_percent)) + geom_point(size = 3) + guides(colour = bar) +
          theme_classic() + labs(title = "scater: total counts vs detected genes", colour = "mito %"))
}
dev.off()
