# # # # # # # # # # # # # # # # # # # # # # # # # 
#
# script for filtering cells 
#
# author: Weizhang
#
# # # # # # # # # # # # # # # # # # # # # # # # # 

Args <- commandArgs()
libCpxInputDir <- Args[6]  # a str, should be seperated by "*"
fragInPeakInput <- Args[7]  # fragInPeakInput output
reportOutput <- Args[8]
figureOutput <- Args[9]
libSizeCutOff <- Args[10]
fragInPeakCutOff <- Args[11]

# extract all libcpx
libcpx <- Sys.glob(file.path(libCpxInputDir, "*.lcpx"))
libcpx.df <- data.frame()
for(i in seq(length(libcpx))){
    libcpx.df[i, 1] <- unlist(strsplit(x = basename(libcpx[i]), split = ".", fixed = TRUE))[1]
    lines <- readLines(con = libcpx[i])
    libcpx.df[i, 2] <- as.numeric(unlist(strsplit(x = lines[8], split = "\t", fixed = TRUE))[3])
}
colnames(libcpx.df) <- c("fileflag", "libcpx")



fragInPeak <- read.table(file = fragInPeakInput)
colnames(fragInPeak) <- c("fileflag", "FragInPeak")
fragInPeak$fileflag <- as.character(x = fragInPeak$fileflag)
for(i in seq(length(fragInPeak$fileflag))){
    fragInPeak$fileflag[i] <- unlist(strsplit(x = basename(fragInPeak$fileflag[i]), split = ".", fixed = TRUE))[1]
}



tmp.fragInPeak <- c()
for(fileflag in libcpx.df$fileflag){
    percent <- fragInPeak[which(fragInPeak$fileflag == fileflag), 2]
    tmp.fragInPeak <- c(tmp.fragInPeak, percent)
}

libcpx.df[["fragInPeak"]] <- tmp.fragInPeak

libcpx.df$libcpx <- log10(as.numeric(libcpx.df$libcpx))
libcpx.df$fragInPeak <- as.numeric(100*libcpx.df$fragInPeak)



# convert
libSizeCutOff <- log10(as.numeric(libSizeCutOff))
fragInPeakCutOff <- 100*as.numeric(fragInPeakCutOff)

options(bitmapType='cairo')
bmp(figureOutput)
plot(libcpx.df$libcpx, libcpx.df$fragInPeak, main = "Cell Filter",
     sub = "select cells locate in topright region",
     xlab = "Library size (log10 reads)", ylab = "Fragments in peaks (%)")
abline(v = libSizeCutOff, col = "black", lwd = 2, lty = 2)
abline(h = fragInPeakCutOff, col = "black", lwd = 2, lty = 2)
dev.off()

idx <- which(libcpx.df$libcpx >= libSizeCutOff & libcpx.df$fragInPeak >= fragInPeakCutOff)

output <- libcpx.df[idx, ]

write.table(x = output, file = reportOutput, quote = FALSE, row.names = FALSE, col.names = FALSE)
