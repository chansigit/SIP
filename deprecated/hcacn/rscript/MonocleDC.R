Args <- commandArgs()
imagedata <- Args[6]
#default 3
num_PCA <- as.numeric(Args[7])
#default 2
cluster_num <- as.numeric(Args[8])
outputpath_now <- Args[9]
print(outputpath_now)
print(Args)
library(monocle)
load(imagedata)
UMIdata <- reduceDimension(UMIdata, max_components=2, num_dim = num_PCA, 
                           reduction_method = 'tSNE', verbose = T) 
UMIdata <- clusterCells(UMIdata, num_clusters=cluster_num, method = "densityPeak")
filename <- as.character(paste(as.character(outputpath_now),"/", "densitypeak_cluster.jpg",sep = ""))
jpeg(filename)
plot_cell_clusters(UMIdata, 1, 2, color_by = "Cluster")
dev.off()
