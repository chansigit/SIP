#!/usr/bin/env Rscript
#Use: Rscript SeruatAnalysis.R --filepath="/path/to/dge/matrix/file"

suppressMessages(library("optparse"))
option_list = list(
  make_option( 
    c("--filepath"),  
    type = "character",
    default = NULL, help = "Dataset file path"
  ),
  
  make_option(
    c("--mincells"),
    type = "numeric",
    default = 1, help = "Include genes with detected expression in at least this many cells"
  ),
  
  make_option(
    c("--mingenes"),
    type = "numeric",
    default = 200, help = "Include cells where at least this many genes are detected"
  ),
  
  make_option(
    c("--nGeneLowerBound"),
    type = "numeric",
    default = 200, help = "Discard cells with nGene < this many detected genes"
  ),
  
  make_option(
    c("--nGeneUpperBound"),
    type = "numeric",
    default = Inf, help = "Discard cells with nGene > this many detected genes"
  ),
  
  make_option(
    c("--mitoLowerBound"),
    type = "numeric",
    default = -Inf, help = "Discard cells with mitochondrial gene ratios < this value"
  ),
  
  make_option(
    c("--mitoUpperBound"),
    type = "numeric",
    default = 0.3, help = "Discard cells with mitochondrial gene ratios > this value"
  ),
  
  make_option(
    c("--resolution"),
    type = "numeric",
    default = 1.03, help = "Clustering parameter: resolution"
  ),
  
  make_option(
    c("--tsneperp"),
    type = "numeric",
    default = 4, help = "t-SNE perplexity"
  ),
  
  make_option(
    c("--DEGeneKept"),
    type = "numeric",
    default = 10, help = "Keep this many DE genes in results"
  ),
  
  make_option(
    c("--DEGeneViz"),
    type = "numeric",
    default = 3, help = "Use this many DE genes for visualization"
  ),
  
  make_option(
    c("--DEminDetectRatio"),
    type = "numeric",
    default = 0.15, help = "DE genes must appear in cells more than this percentage"
  ),
  
  make_option(
    c("--DEpvalcut"),
    type = "numeric",
    default = 0.05, help = "DE gene detection P-value"
  )
);

opt_parser = OptionParser(option_list=option_list);
params= parse_args(opt_parser);

if (is.null(params$filepath)){
  print_help(opt_parser)
  stop("File path must be supplied (input file).\n", call.=FALSE)
}


########################################################################################
########################################################################################
# Data Loading Params

param.min.genes <- 200

# QC and Feature Engineering Params

# Clustering and Visualization Params

# Finding markers

suppressMessages(library(Seurat))
expressionMatrix <- read.table(params$filepath)
# ----------------------------------------------------------------------------------------
#                                     Loading Data
##  https://learn.gencore.bio.nyu.edu/single-cell-rnaseq/seurat-part-1-loading-the-data/
##  min.cells := Include genes with detected expression in at least this many cells
##  min.genes := Include cells where at least this many genes are detected
# ----------------------------------------------------------------------------------------
cellhub <- CreateSeuratObject(
  project   = "SIP_Seurat_Analysis",
  raw.data  =  expressionMatrix,
  min.cells =  params$mincells,
  min.genes =  params$mingenes
)

show.compression=FALSE
if (show.compression){
  dense.size <- object.size(x = as.matrix(x = cellhub@data))
  sparse.size <- object.size(x = cellhub@data) 
  print(paste("Dense size=",dense.size/1024.0/1024.0,"MB" ))
  print(paste("Spare size=",sparse.size/1024.0/1024.0,"MB"))
  print(paste("Saving=",as.numeric(dense.size/sparse.size)," x"))
}




# ----------------------------------------------------------------------------------------
#                               Mitochondrial Gene QC
## So mito gene rate is important QC metric，it can't be too high. 
## We obtain the mitochondrial gene list in this way:
##     suppressMessages(library(EnsDb.Hsapiens.v86))
##     gns <- genes(EnsDb.Hsapiens.v86, filter = ~ seq_name == "MT")
##     mito.gene.db <- rownames(as.data.frame(gns))
# ----------------------------------------------------------------------------------------
## construct mitochondrial gene list
mito.gene.db<-c(
  'ENSG00000210049', 'ENSG00000211459', 'ENSG00000210077', 'ENSG00000210082', 
  'ENSG00000209082', 'ENSG00000198888', 'ENSG00000210100', 'ENSG00000210107',
  'ENSG00000210112', 'ENSG00000198763', 'ENSG00000210117', 'ENSG00000210127',
  'ENSG00000210135', 'ENSG00000210140', 'ENSG00000210144', 'ENSG00000198804',
  'ENSG00000210151', 'ENSG00000210154', 'ENSG00000198712', 'ENSG00000210156',
  'ENSG00000228253', 'ENSG00000198899', 'ENSG00000198938', 'ENSG00000210164',
  'ENSG00000198840', 'ENSG00000210174', 'ENSG00000212907', 'ENSG00000198886',
  'ENSG00000210176', 'ENSG00000210184', 'ENSG00000210191', 'ENSG00000198786',
  'ENSG00000198695', 'ENSG00000210194', 'ENSG00000198727', 'ENSG00000210195', 
  'ENSG00000210196')
mito.genes <- c()

## The gene ID we use in the expression matrix is the with-version ENSG ID 
## Here we map the no-version IDs to the with-version IDs.
for (mit_ensg in mito.gene.db){ # for each gene id in annotated mito gene db
  regex=paste0("^",mit_ensg)  # find the corresponding genes
  gene<-grep(pattern=regex, x=rownames(x=cellhub@data), value=TRUE) 
  mito.genes<- c(mito.genes, gene)
}
#print(mito.genes)


# Calculate the percentage of mitochondrial genes here and store it in Seurat's metadata
suppressMessages(library(Matrix))
percent.mito <- Matrix::colSums(cellhub@raw.data[mito.genes, ])/ # NOTE: load the Matrix package 
  Matrix::colSums(cellhub@raw.data)
cellhub <- AddMetaData(object = cellhub, 
                       metadata = percent.mito, 
                       col.name = "MitoGeneRatio")
#percent.mito



# ----------------------------------------------------------------------------------------
#                               Expression-level QC
# The number of genes and UMIs (nGene and nUMI) are automatically calculated by Seurat.
# Cells having a low number of detected genes might represent low quality cells.
# nGene is the number of genes detected (a gene is called detected if its expression level
# passes the min.cells and min.genes filtering in the data loading process)
# For non-UMI data, nUMI represents the sum of the non-normalized values within a cell 
# Here we plot the violin plots for nGene, nUMI, and Mitochondrial percentages.
# ----------------------------------------------------------------------------------------
pdf("CellQC.pdf")
par(mfrow = c(2, 2))
par(mar = c(4, 4, 1, 8) + 0.1)
VlnPlot(object = cellhub, features.plot = c("nGene", "nUMI", "MitoGeneRatio"), nCol = 2)
dev.off()


# ----------------------------------------------------------------------------------------
#                                   Do Cell QC
## Too low nGene or nUMI indicates lysing cells, while too high indicates doubletsression, 
## High level of mitochondrial gene expression indicates: 1) poor sample quality, apoptotic
## or lysing cells; 2) special samples, some tumor biopsies have such characteristic.
## see ref. https://bit.ly/2DAj3gC
# ----------------------------------------------------------------------------------------
cellhub <- FilterCells(cellhub,
                       subset.names   = c("nGene", "MitoGeneRatio"), 
                       low.thresholds = c(params$nGeneLowerBound,params$mitoLowerBound), 
                       high.thresholds= c(params$nGeneUpperBound,params$mitoUpperBound)
)


# More information about seurat slots can be found here:
# https://github.com/satijalab/seurat/wiki/Seurat#slots
#print(slotNames(cellhub))
#cellhub


# ----------------------------------------------------------------------------------------
#                                     Normalization
# By default, Seurat implements a global-scaling normalization method “LogNormalize” 
# that normalizes the gene expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
# ----------------------------------------------------------------------------------------
cellhub <- NormalizeData(
  cellhub, display.progress=F,
  normalization.method="LogNormalize", 
  scale.factor=10000
)


# ----------------------------------------------------------------------------------------
#                                  Feature Selection 
## 1) Calculate average expresion and dispersion of each gene
##           gene1.avgexp   gene2.avgexp  ...  gene2w.avgexp
##           gene1.var      gene2.var     ...  gene2w.var
## 2) Divide genes into 20 bins based on average expression
##      [gene1, gene2,..., gene1000]        [gene1001, gene1002,..., gene2000]
##      [gene2001, gene2002,..., gene3000]  [gene3001, gene3002,..., gene4000]
##                    ...                   [gene19991,gene19992,...,gene20000]
## 3) Calculate  z-scores for dispersion within each bin.
##        z-score([gene1.var,     gene2.var,...,     gene1000.var])    
##                                ...
##        z-score([gene19991.var, gene19992.var,..., gene20000.var])    
# ----------------------------------------------------------------------------------------
pdf("MeanDispersionPlot.pdf")
cellhub <- FindVariableGenes(object = cellhub,  
                             mean.function = ExpMean, dispersion.function = LogVMR, 
                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,
                             do.plot=TRUE, plot.both=FALSE, do.text=FALSE, contour.lwd=0.5, contour.col="red" ,contour.lty = 3, cex.use=0.1                       
)
dev.off()
print(paste0("VarGene#=",length(x=cellhub@var.genes)))
#par(mfrow = c(1, 2))
#GenePlot(object = cellhub, gene1 = "nUMI", gene2 = "MitoGeneRatio")
#GenePlot(object = cellhub, gene1 = "nUMI", gene2 = "nGene")


# ----------------------------------------------------------------------------------------
#                          Scale Batch Effect & Confounding Factors 
## Scales and centers genes in the dataset
## Regress out cell-cell variation in gene expression driven by batch (if any)
## Batch effect removal can be based on 1) nUMI 2) MitoGeneRatio 3) cell
## alignment rate
## Other confounding factors can be removed by learned "cell-cycle" score.
# ----------------------------------------------------------------------------------------
cellhub <- ScaleData(object = cellhub,  do.cpp=TRUE, do.par=TRUE,
                     vars.to.regress = c("nUMI", "MitoGeneRatio"))

# ----------------------------------------------------------------------------------------
#                              PCA dimensionality reduction
# ----------------------------------------------------------------------------------------
cellhub <- RunPCA(object   = cellhub, pc.genes = cellhub@var.genes, 
                  do.print = FALSE, pcs.print = 1:12, genes.print = 10)

pdf("VizPCA.pdf")
par(mar = c(4, 8, 1, 2) + 0.1)
VizPCA(object = cellhub, pcs.use =1:6, nCol=2)
VizPCA(object = cellhub, pcs.use =7:12, nCol=2)
dev.off()
# Evaluating PCA performance
# cellGroupObj <- JackStraw(object = cellGroupObj, num.replicate =    100, display.progress = FALSE)
# JackStrawPlot(object=cellGroupObj, PCs=1:12)

# ----------------------------------------------------------------------------------------
#                              Sub-population identification
##  Seurat used a graph-based clustering to indentify sub-groups. The method build a KNN
##  graph based on Euclidean distance in PCA space and find clique on the graph.  
##  **Resolution** is an important parameter in this population identification step.
##  Value of the resolution parameter:
##      if you want to obtain a larger  number of communities.
##      use a value > 1.0 
##      if you want to obtain a smaller number of communities.
##      use a value < 1.0 
# ----------------------------------------------------------------------------------------
cellhub <- FindClusters(object=cellhub, save.SNN=TRUE,
                        reduction.type="pca", dims.use=1:12,
                        plot.SNN= F, force.recalc =T, print.output=F,
                        resolution = params$resolution
)

groupinfo <- as.data.frame(cellhub@ident)
#groupinfo
write.table(groupinfo, "GroupInfo.tsv", sep="\t", col.names=NA)



# ----------------------------------------------------------------------------------------
#                               t-SNE visualization
## t-SNE is a non-linear method to visualize high-dimensional data
## The **perplexity** parameter is very important.
## For more info about t-SNE https://distill.pub/2016/misread-tsne/
# ----------------------------------------------------------------------------------------
cellhub <- RunTSNE(object = cellhub, dims.use = 1:12, do.fast = TRUE,
                   perplexity= params$tsneperp)

pdf("tSNE.pdf")
TSNEPlot(object=cellhub)#
dev.off()

save(cellhub, file="SeuratObjectSerialized.Robj")

# ----------------------------------------------------------------------------------------
#                           Find Marker Genes
## FindAllMarkers reports DE genes for every cluster compared to all
## remaining cells, (report the positive ones only). 
## **min.pct** argument requires a gene to be detected at a minimum percentage in either
##   of the two groups of cells. 
## **return.thres** parameter sets the reporting cut-off 
# ----------------------------------------------------------------------------------------
markerGenes <- FindAllMarkers(
  object = cellhub, 
  test.use = "bimod",
  only.pos = TRUE, 
  min.pct  = params$DEminDetectRatio, 
  return.thresh = params$DEpvalcut
)



suppressMessages(library(dplyr))
# use of %>% requires dplyr
topDEGenes <- markerGenes %>% group_by(cluster) %>% top_n(params$DEGeneKept, avg_logFC)  
vizDEGenes <- markerGenes %>% group_by(cluster) %>% top_n(params$DEGeneViz,  avg_logFC)
vizDEGenes <- as.vector(vizDEGenes[["gene"]])

write.table(topDEGenes, "DE_Genes.tsv", sep="\t", col.names=NA)

pdf("DE_Genes_ViolinPlot.pdf")
VlnPlot(
  object = cellhub, 
  features.plot = vizDEGenes, 
  size.title.use=10, size.x.use=10
)
dev.off()

pdf("DE_Genes_FeaturePlot.pdf")
FeaturePlot(
  object = cellhub, 
  features.plot = vizDEGenes, 
  cols.use = c("grey", "blue"), 
  reduction.use = "tsne", dark.theme=FALSE,
  vector.friend=TRUE, pt.size=10
)
dev.off()

pdf("DE_Genes_Heatmap.pdf")
DoHeatmap(object = cellhub, 
          genes.use = vizDEGenes, 
          slim.col.label = TRUE, remove.key = TRUE
)
dev.off()
