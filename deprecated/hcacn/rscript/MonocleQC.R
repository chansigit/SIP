Args <- commandArgs()
matrixdata <- Args[6]
#default 0.1
min_expression<- as.numeric(Args[7])
#default 1
num_cells_expressed_threshold<- as.numeric(Args[8])
#default 1e6
TotalmRNAs <- as.numeric(Args[9])
#default 0.1
mean_expression_threshold <- as.numeric(Args[10])
outputpath <- Args[11]
library(monocle)
#Read in data
#rawdata <- read.table("out_gene_exon_tagged.txt", row.names = 1,header = T ,stringsAsFactors = F)
rawdata <- read.table(matrixdata, row.names = 1,header = T ,stringsAsFactors = F)
#Construct pd and fd. This step is meaningless, just construct required inputs for newCellDataSet
sample_sheet <- data.frame(groups = rep(1:2, times=c(1,ncol(rawdata)-1)), row.names = colnames(rawdata))
gene_ann <- data.frame(gene_short_name = row.names(rawdata), row.names = row.names(rawdata))
pd <- new("AnnotatedDataFrame",data=sample_sheet)
fd <- new("AnnotatedDataFrame",data=gene_ann)
#Build CellDataSet object
UMIdata <- newCellDataSet(as.matrix(rawdata),phenoData = pd,featureData =fd,
                          expressionFamily = negbinomial.size())
#Estimate size factors and dispersions 
UMIdata <- estimateSizeFactors(UMIdata)
UMIdata <- estimateDispersions(UMIdata)
#Filtering low-quality cells,"min_expr" is the expression threshold
UMIdata <- detectGenes(UMIdata, min_expr = min_expression)
#print(head(fData(UMIdata)))
#"num_cells_expressed" means at least "num_cells_expressed" cells express these genes
expressed_genes <- row.names(subset(fData(UMIdata),num_cells_expressed >= num_cells_expressed_threshold))
#print(head(pData(UMIdata)))
pData(UMIdata)$Total_mRNAs <- Matrix::colSums(exprs(UMIdata))
UMIdata <- UMIdata[,pData(UMIdata)$Total_mRNAs < TotalmRNAs]
upper_bound <- 10^(mean(log10(pData(UMIdata)$Total_mRNAs)) +
                     2*sd(log10(pData(UMIdata)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(UMIdata)$Total_mRNAs)) -
                     2*sd(log10(pData(UMIdata)$Total_mRNAs)))
jpeg(filename = as.character(paste(as.character(outputpath),"/", "density_Total_mRNAs.jpg",sep = "")))
qplot(Total_mRNAs, data = pData(UMIdata), geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)
dev.off()
disp_table <- dispersionTable(UMIdata)
#head(disp_table)
#Select >="mean_expression" genes
unsup_clustering_genes <- subset(disp_table, mean_expression >= mean_expression_threshold)
UMIdata <- setOrderingFilter(UMIdata, unsup_clustering_genes$gene_id)
#Plot the curve of mean_expression and disperssion_emperical
#The red curve shows the mean-variance model learning by estimateDispersions().
jpeg(filename = as.character(paste(as.character(outputpath),"/", "meanexpression_disersionemprical.jpg",sep = "")))
plot_ordering_genes(UMIdata)
dev.off()
jpeg(filename = as.character(paste(as.character(outputpath),"/", "PCvariance.jpg",sep = "")))
plot_pc_variance_explained(UMIdata, return_all = F,max_components = 15) 
dev.off()
save.image(as.character(paste(as.character(outputpath),"/", "MonocleQCimage.Rdata",sep = "")))



