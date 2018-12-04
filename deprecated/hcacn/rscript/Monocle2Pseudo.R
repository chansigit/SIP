Args <- commandArgs()
print(Args)
cds_file <- Args[6]

result_folder <- Args[7]

# default "tSNE"
# reduction_method<-
# #default 3

# num_reduce_dimension <- as.numeric(Args[7])
# #default 2
# num_max_components <- as.numeric(Args[8])
# outputpath_now <- Args[9]

library(monocle)
load(cds_file)


fData(cds)$num_cells_expressed <-apply(exprs(cds)==0,1,sum)
expressed_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))

cds <- reduceDimension(cds,
                      max_components = 2,
                      norm_method = 'log',
                      num_dim = 3,
                      reduction_method = 'tSNE',
                      verbose = T)
cds <-clusterCells(cds,
                 rho_threshold = 2,
                 delta_threshold = 4,
                 skip_rho_sigma = T,
                 verbose = F)


clustering_DEG_genes <- differentialGeneTest(cds[expressed_genes,],
									          fullModelFormulaStr = '~Cluster',
									          cores = 24)
ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
cds <-setOrderingFilter(cds,ordering_genes = ordering_genes)
cds <-reduceDimension(cds, method = 'DDRTree')
cds <-orderCells(cds)

print("begin Plotting")
jpeg(filename = paste(result_folder, "ordering_genes.jpg",sep = "/"))
plot_ordering_genes(cds)
dev.off()

jpeg(filename = paste(result_folder,"cell_clusters.jpg",sep = "/"))
plot_cell_clusters(cds, color_by = 'as.factor(Cluster)')
dev.off()

jpeg(filename = paste(result_folder, "cell_trajectory_by_groups.jpg",sep = "/"))
plot_cell_trajectory(cds,color_by = 'groups')
dev.off()

jpeg(filename = paste(result_folder, "cell_trajectory_by_Pseudotime.jpg",sep = "/"))
plot_cell_trajectory(cds,color_by = 'Pseudotime')
dev.off()

jpeg(filename = paste(result_folder, "cell_trajectory_by_Cluster.jpg",sep = "/"))
plot_cell_trajectory(cds,color_by = 'Cluster')
dev.off()

jpeg(filename = paste(result_folder, "cell_trajectory_by_State.jpg",sep = "/"))
plot_cell_trajectory(cds,color_by = 'State')
dev.off()


# save cds object
save(cds,file=paste(result_folder,'cds_Pseudo.Rdata',sep = "/"))


