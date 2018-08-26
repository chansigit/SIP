library("SingleCellExperiment")
Args <- commandArgs()
print(Args)

rawdata_file <- Args[6]
ann_file <-Args[7]
# LOG,ORIGIN
rawdata_format <- Args[8]
RData_file <- Args[9]



matrix<-read.table(rawdata_file,sep = ",",row.names=1,header = TRUE,)

  ann <- tryCatch({
    read.csv(ann_file)
  }, error = function(e) {
    data.frame(row.names = colnames(matrix))
  })
show("Annotation File dimension:")
show(dim(ann))
show("Matrix File dimension:")
show(dim(matrix))

# Generate sce 类
if(rawdata_format == "ORIGIN")
    {
        sce <- SingleCellExperiment(
            assays = list(
                counts = as.matrix(matrix),
                logcounts = log2(as.matrix(matrix) + 1)
                ), 
                colData = ann
            )
        # define feature names in feature_symbol column
        rowData(sce)$feature_symbol <- rownames(sce)

        ### preprocess!
        # remove features with duplicated names
        sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
        # define spike-ins
        isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)


    }else if (rawdata_format =="LOG")    {
        show("Have not implementation till now...")
    }else     {
        show("Wrong raw data Format!")
    }

save(sce,file=RData_file)

