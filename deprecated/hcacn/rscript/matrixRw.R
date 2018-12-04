Args <- commandArgs()
matrixdata <- Args[6]
outputpath <- Args[7]
rawmatrix <- read.table(file = matrixdata, header = T, row.names = 1)
write.table(rawmatrix, file = as.character(paste(as.character(outputpath),"/", "processedMatrix.txt",sep = "")))