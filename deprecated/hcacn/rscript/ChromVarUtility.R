# # # # # # # # # # # # # # # # # # # # # # # # # 
#
# script for running chromVAR
#
# author: Weizhang
#
# # # # # # # # # # # # # # # # # # # # # # # # # 

library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(Matrix)
library(pheatmap)
library(BiocParallel)
set.seed(2018)


Args <- commandArgs()
bamDir <- Args[6]  # bamfile dictionary
peakfile <- Args[7]  # peakfile path
thread <- Args[8]  # number of core
genomeidx <- Args[9]  # specify a genome
var.tiff <- Args[10]  # pdf file path
clu.tiff <- Args[11]  # pdf file path
devmatrix_path <- Args[12]  # devmatrix file path
variability_path <- Args[13]  # variability matrix file path
clu_matrix_path <- Args[14]  # clustring matrix file path

if(genomeidx == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
    genome <- BSgenome.Hsapiens.UCSC.hg19
    motifs <- getJasparMotifs(species = "Homo sapiens")
    # print("BSgenome.Hsapiens.UCSC.hg19")
}else if(genomeidx == "hg38"){
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- BSgenome.Hsapiens.UCSC.hg38
    motifs <- getJasparMotifs(species = "Homo sapiens")
    # print("BSgenome.Hsapiens.UCSC.hg38")
}else if(genomeidx == "mm9"){
    library(BSgenome.Mmusculus.UCSC.mm9)
    genome <- BSgenome.Mmusculus.UCSC.mm9
    motifs <- getJasparMotifs(species = "Mus musculus")
    # print("BSgenome.Mmusculus.UCSC.mm9")
}else if(genomeidx == "mm10"){
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome <- BSgenome.Mmusculus.UCSC.mm10
    motifs <- getJasparMotifs(species = "Mus musculus")
    # print("BSgenome.Mmusculus.UCSC.mm10")
}

# extract all bam files
bamfile <- Sys.glob(file.path(bamDir, "*.bam"))


# peak file
peaks <- getPeaks(peakfile, sort_peaks = TRUE)

# set core
register(SnowParam(workers = 4, type = "SOCK"))

# Counts
colData = S4Vectors::DataFrame(celltype = as.character(seq(length(bamfile))))
example_counts <- getCounts(bamfile, peaks, paired =  TRUE, by_rg = FALSE,
                            format = "bam", colData = colData)

example_counts <- addGCBias(example_counts,
                            genome = genome)

counts_filtered <- filterSamples(example_counts, min_depth = 1500, min_in_peaks = 0.15, shiny = FALSE)
counts_filtered <- filterPeaks(counts_filtered)

motif_ix <- matchMotifs(motifs, counts_filtered, genome = genome)

# computing deviations and save deviation matrix
dev <- computeDeviations(object = counts_filtered, annotations = motif_ix)
devmatrix <- as.data.frame(deviations(dev))
write.table(x = devmatrix, file = devmatrix_path, quote = FALSE, 
            sep = "\t", row.names = TRUE, col.names = TRUE)
# computing variability and save variability
variability <- computeVariability(dev)
variability_matrix <- as.data.frame(cbind(variability$variability, variability$p_value, variability$p_value_adj))

rownames(variability_matrix) <- as.character(variability$name)
colnames(variability_matrix) <- c("variability", "p_value", "p_value_adj")

write.table(x = variability_matrix, file = variability_path, quote = FALSE, 
            sep = "\t", row.names = TRUE, col.names = TRUE)

# linux
options(bitmapType='cairo')

# plot variability
bmp(var.tiff)
plotVariability(variability, use_plotly = FALSE)
dev.off()

# clustering
sample_cor <- getSampleCorrelation(dev)
pheatmap(as.dist(sample_cor),
         annotation_row = colData(dev),
         clustering_distance_rows = as.dist(1-sample_cor),
         clustering_distance_cols = as.dist(1-sample_cor),
         filename = clu.tiff)

# save clu matrix
sample_cor_matrix <- as.data.frame(sample_cor)
write.table(x = sample_cor_matrix, file = clu_matrix_path, quote = FALSE, 
            sep = "\t", row.names = TRUE, col.names = TRUE)


