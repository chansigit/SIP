Args <- commandArgs()
print(Args)

# matrix <- Args[6]
# barcodes <- Args[7]
# genes <- Args[8]
#projectname <- Args[3]
input <- Args[6]
output <- Args[7]
if (!dir.exists(output))
{
	dir.create(output)
}
setwd(output)

library(Seurat)
library(dplyr)
library(Matrix)

# celldata = Matrix::readMM(file = matrix)
# cellname = read.table(file = barcodes, header = FALSE, colClasses = "character")[[1]]
# genename = read.table(file = genes, header = FALSE, colClasses = "character")[[1]]
# rownames(celldata) = genename
# colnames(celldata) = cellname

#load data

celldata = Read10X(data.dir = input)
test = CreateSeuratObject(raw.data = celldata,project = "Frankie so tired")

#processing
#mito.genes <- grep(pattern = "^MT-", x = rownames(x = test@data), value = TRUE)
# percent.mito <- Matrix::colSums(test@raw.data[mito.genes, ])/Matrix::colSums(test@raw.data)
# test <- AddMetaData(object = test, metadata = percent.mito, col.name = "percent.mito")

##############################################################################################
#file.create("violinplot.jpeg")
jpeg(file="violinplot.jpeg")
VlnPlot(object = test, features.plot = c("nGene", "nUMI"), nCol = 2)
dev.off()

##############################################################################################
#file.create("geneplot.jpeg")
jpeg(file="geneplot.jpeg")
GenePlot(object = test, gene1 = "nUMI", gene2 = "nGene")
dev.off()

#
test <- FilterCells(object = test, subset.names = c("nGene","nUMI"),
                  low.thresholds = c(-Inf, -Inf), high.thresholds = c(5600,40000))


test <- NormalizeData(object = test, normalization.method = "LogNormalize",scale.factor = 10000)

##############################################################################################
#file.create("variableGenes.jpeg")
jpeg(filename ="variableGenes.jpeg")
test <- FindVariableGenes(object = test, mean.function = ExpMean, dispersion.function = LogVMR,
                            x.low.cutoff = 0.05, x.high.cutoff = 3, y.cutoff = 0.5)
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
#
# # FeaturePlot(test,features.plot = "S100B",cols.use = c("white","red"))
# # cluster3.markers <- FindMarkers(object = test, ident.1 = 1, logfc.threshold = 0.8)
# # print(x = head(x = cluster4.markers, n = 5))

 