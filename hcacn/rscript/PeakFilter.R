# please modify before using!!!!!!!!!!
# script for discard blacklist region
# run in the current dictionary
# need packages: rtracklayer

Args <- commandArgs()
User_bed_file <- Args[6]  # this is user's peak summit file from macs2
blacklist_bed_file <- Args[7]  # this is blacklist file from UCSC or user
output_bed_file <- Args[8]  # this is output file
overlap_rate <- as.numeric(Args[9])  # this is overlap rate for discard peak
extendRange <- as.numeric(Args[10])  # this is extend range up/downstream
top <- as.numeric(Args[11])  # extract top number peaks using parameter top, if is 0, means all peaks

library(rtracklayer)
library(GenomicRanges)

gr_a <- rtracklayer::import(con = blacklist_bed_file, format = "bed")

gr_b <- rtracklayer::import(con = User_bed_file, format = "bed")
start(gr_b) <- start(gr_b) - extendRange
end(gr_b) <- end(gr_b) + extendRange

o <- GenomicRanges::findOverlaps(query = gr_a, subject = gr_b, ignore.strand = TRUE)
o_1 <- gr_a[S4Vectors::queryHits(o)]
o_2 <- gr_b[S4Vectors::subjectHits(o)]
peak_intersect <- IRanges::pintersect(o_1, o_2)
discard_flag <- (IRanges::width(peak_intersect)/pmin(IRanges::width(o_1),IRanges::width(o_2)) >= overlap_rate)
discard_reads <- IRanges::unique(o_2[discard_flag])
remain_reads <- gr_b[!gr_b %in% discard_reads]
if(top < 0){
    stop("Parameter top must >= 0!!!")
}else if(top == 0){
    rtracklayer::export(object = remain_reads, con = output_bed_file, format = "BED")
}else{
    gr_a <- sort(remain_reads, by = ~ score, decreasing = TRUE)
    gr_a <- sort(gr_a[1:top])
    print("Extracting peaks......")
    rtracklayer::export(object = gr_a, con = output_bed_file, format = "BED")
}



