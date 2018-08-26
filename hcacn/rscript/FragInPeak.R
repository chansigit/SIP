# # # # # # # # # # # # # # # # # # # # # # # # # 
#
# script for calculating percentage of fragments in peaks
#
# author: Weizhang
#
# # # # # # # # # # # # # # # # # # # # # # # # # 

library(rtracklayer)

Args <- commandArgs()
fragDir <- Args[6]  # fragments bed dir
peakfile <- Args[7]  # peakfile path
output_path <- Args[8]  # output file path

# extract all bed files
fragfile <- Sys.glob(file.path(fragDir, "*.bed"))
fragfile <- sort(fragfile)

gr_peak <- import(con = peakfile, format = "BED")

tbOutput <- data.frame()

for(i in seq(length(fragfile))){
    gr_frag <- import(con = fragfile[i], format = "BED")
    frag.num <- length(gr_frag)
    o <- GenomicRanges::findOverlaps(query = gr_frag, subject = gr_peak, ignore.strand = TRUE)
    fraginpeak.num <- length(unique(S4Vectors::queryHits(o)))
    fraginpeak.percentage <- fraginpeak.num/frag.num
    tbOutput[i, 1] <- fragfile[i]
    tbOutput[i, 2] <- fraginpeak.percentage
}

write.table(x = tbOutput, file = output_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


