# Created by: Fackie
# Created on: 2018/3/15
Args <- commandArgs()
print(Args)

input <- Args[6]
output <- Args[7]
xlow <- Args[8]
xhigh <- Args[9]
y <- Args[10]

if (!dir.exists(output))
{
  dir.create(output)
}
#setwd(output)
library(Seurat)
library(dplyr)
library(Matrix)
load(input)
setwd(output)
#############################################################################################
# humi <- quantile(test@meta.data$nUMI,0.95)
# hgene <- quantile(test@meta.data$nGene,0.95)
# lumi <- quantile(test@meta.data$nUMI,0.05)
# lgene <- quantile(test@meta.data$nGene,0.05)
#
# test <- FilterCells(object = test, subset.names = c("nGene","nUMI"), low.thresholds = c(lgene, lumi), high.thresholds = c(hgene,humi))
#
# test <- NormalizeData(object = test, normalization.method = "LogNormalize",scale.factor = 10000)

if(xlow == 'xlow' & xhigh == 'xhigh' && y == 'y')
{
  xhigh <- quantile(test@hvg.info$gene.mean,0.9)
  xlow <- quantile(test@hvg.info$gene.mean,0.1)
  y <- quantile(test@hvg.info$gene.dispersion.scaled,0.9)
}
##############################################################################################
#file.create("variableGenes.jpeg")
jpeg(filename ="variableGenes.jpeg")
test <- FindVariableGenes(object = test, mean.function = ExpMean, dispersion.function = LogVMR,
                          x.low.cutoff = xlow, x.high.cutoff = xhigh, y.cutoff = y)
dev.off()


test <- ScaleData(object = test,vars.to.regress = "nUMI")
test <- RunPCA(object = test, pc.genes = test@var.genes)

##############################################################################################
#file.create("Elbowplot.jpeg")
jpeg(file="Elbowplot.jpeg")
PCElbowPlot(object = test,num.pc = 30)
dev.off()

test <- FindClusters(object = test, reduction.type = "pca", print.output = 1)
test <- RunTSNE(object = test,dims.use = 1:10, perplexity =10, do.fast = TRUE)

##############################################################################################
#file.create("TSNEplot.jpeg")
jpeg(file="TSNEplot.jpeg")
TSNEPlot(object = test)
dev.off()
