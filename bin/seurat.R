#!/usr/bin/env Rscript
# SIP downstream with Seurat: the paper's 9-step workflow (cluster finding + DEG analysis).
# Usage: seurat.R matrix=<tsv> genes=<tsv> key=value ...
suppressMessages({ library(Seurat); library(ggplot2) })
a <- commandArgs(TRUE)
p <- setNames(as.list(sub("^[^=]*=", "", a)), sub("=.*", "", a))
num <- function(k) as.numeric(p[[k]])
set.seed(1)
figs <- list()   # every figure, also handed to the summary report (seurat_plots.rds)
out <- function(file, plots, width = 7, height = 7) {
    figs[names(plots)] <<- plots
    pdf(file, width = width, height = height); for (x in plots) print(x); dev.off()
}
bar <- guide_colourbar(display = "rectangles")   # vector colour bars instead of an embedded gradient image
# DoHeatmap still embeds two bitmaps: the cluster bar (annotation_raster) and the colour bar set inside the
# fill scale. Redraw the cluster bar as rectangles and swap the guide, so the figure is pure vector.
vector_heatmap <- function(p) {
    for (k in which(vapply(p$layers, function(l) inherits(l$geom, "GeomRasterAnn"), TRUE))) {
        a <- p$layers[[k]]$geom_params; r <- as.matrix(a$raster); n <- ncol(r); w <- (a$xmax - a$xmin) / n
        p$layers[[k]] <- annotate("rect", xmin = a$xmin + (seq_len(n) - 1) * w, xmax = a$xmin + seq_len(n) * w,
                                  ymin = a$ymin, ymax = a$ymax, fill = r[1, ], colour = NA)
    }
    fill <- p$scales$get_scales("fill"); fill$guide <- bar
    p
}

m <- as.matrix(read.delim(p$matrix, row.names = 1, check.names = FALSE))
g <- read.delim(p$genes, header = FALSE, row.names = 1)
rownames(m) <- make.unique(ifelse(rownames(m) %in% rownames(g), g[rownames(m), 1], rownames(m)))

# 1) initialize a Seurat object
so <- CreateSeuratObject(m, min.cells = num("min_cells"), min.features = num("min_features"))
mt <- grep("^MT-", rownames(so), ignore.case = TRUE, value = TRUE)
so$percent.mt <- 0
if (length(mt)) so[["percent.mt"]] <- PercentageFeatureSet(so, features = mt)
out("cell_qc.pdf", list("Seurat: cell QC" = VlnPlot(so, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)))

# 2) filter cells on built-in metrics
so <- subset(so, percent.mt <= num("max_mt"))
# 3) normalize
so <- NormalizeData(so)
# 4) select variable genes
so <- FindVariableFeatures(so)
out("variable_genes.pdf", list("Seurat: variable genes" = VariableFeaturePlot(so)))
# 5) remove confounding factors by regression
rv <- strsplit(p$regress_vars, ",")[[1]]
so <- ScaleData(so, features = rownames(so), vars.to.regress = if (length(rv)) rv)
# 6) linear dimensionality reduction
npc <- min(num("n_pcs"), ncol(so) - 1)
so <- RunPCA(so, npcs = npc, verbose = FALSE)
out("pca.pdf", list("Seurat: PCA loadings" = VizDimLoadings(so, dims = 1:min(6, npc), ncol = 2, nfeatures = 15),
                   "Seurat: PCA elbow plot" = ElbowPlot(so, ndims = npc)), 10, 12)
# 7) t-SNE (perplexity capped by the cell number)
so <- RunTSNE(so, dims = 1:min(num("tsne_dims"), npc),
              perplexity = min(num("tsne_perplexity"), floor((ncol(so) - 1) / 3)))
# 8) PCA-based clustering
cd <- 1:min(num("cluster_dims"), npc)
so <- FindNeighbors(so, dims = cd, k.param = min(20, ncol(so) - 1))
so <- FindClusters(so, resolution = num("resolution"))
out("tsne.pdf", list("Seurat: t-SNE, coloured by cluster" = DimPlot(so, reduction = "tsne", pt.size = 3)))
write.table(data.frame(cell = colnames(so), cluster = Idents(so)), "clusters.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 9) differentially expressed genes
if (nlevels(Idents(so)) > 1) {
    mk <- FindAllMarkers(so, only.pos = TRUE, min.pct = num("de_min_pct"))
    write.table(mk, "markers.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    n <- ceiling(num("heatmap_genes") / nlevels(Idents(so)))
    top <- unlist(lapply(split(mk, mk$cluster), function(d) head(d$gene[order(-d$avg_log2FC)], n)))
    top <- head(unique(top), num("heatmap_genes"))
    hm <- vector_heatmap(DoHeatmap(so, features = top, raster = FALSE))
    out("de_genes.pdf", list(
        "Seurat: marker gene heatmap" = hm,
        "Seurat: top marker genes" = VlnPlot(so, head(top, 6), ncol = 3),
        "Seurat: top marker genes on t-SNE" = FeaturePlot(so, head(top, 6), reduction = "tsne", ncol = 3) & guides(colour = bar)), 10, 8)
}
saveRDS(so, "seurat.rds")
saveRDS(figs, "seurat_plots.rds")
