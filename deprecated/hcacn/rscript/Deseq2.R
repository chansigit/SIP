 Args <- commandArgs()
 inputmatrix <- Args[6]
 inputcondition <- Args[7]
 outputpath <- Args[8]
 p <- as.numeric(Args[9])
 lgFDl <- as.numeric(Args[10])
 lgFDu <- as.numeric(Args[11]) 
 library(DESeq2)
 library(pheatmap)
 expressionmatrix <- read.table(file = inputmatrix, header = T, row.names = 1)
 condition <- read.csv(file = inputcondition,stringsAsFactors = F, header = T)
 condition$condition<-factor(condition$condition, ordered = F)
 dds <- DESeqDataSetFromMatrix(expressionmatrix, colData = condition, design= ~ condition)
 dds <- DESeq(dds)  
 res <- results(dds)
 res <- res[order(res$padj),]
 sink(file = as.character(paste(as.character(outputpath),"/", "Summary.txt",sep = "")))
 res
 sink()
 diff_gene_deseq2 <- subset(res, padj < p & (log2FoldChange > lgFDu | log2FoldChange < lgFDl))
 diff_gene_deseq2 <- as.data.frame(diff_gene_deseq2)
 jpeg(filename = as.character(paste(as.character(outputpath),"/", "MAplot.jpg",sep = "")))
 plotMA(dds, main="DESeq2_MAplot", ylim=c(min(as.vector(res$log2FoldChange), na.rm = T),max(as.vector(res$log2FoldChange), na.rm = T)))
 dev.off()
 select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
 # defaults to log2(x+1)
 nt <- normTransform(dds) 
 log2.norm.counts <- assay(nt)[select,]
 df <- as.data.frame(colData(dds)[,c("name","condition")])
 png(filename = as.character(paste(as.character(outputpath),"/", "DEheatmap.png",sep = "")), width=800, height = 1000)
 pheatmap(log2.norm.counts,color = colorRampPalette(c("navy", "white", "firebrick3"))(50), cluster_rows=TRUE, show_rownames=T,
          cluster_cols=TRUE, annotation_col=df, scale = "row")
 dev.off()
 write.csv(diff_gene_deseq2, file = as.character(paste(as.character(outputpath),"/", "Deseq2Result.csv",sep = "")))
 