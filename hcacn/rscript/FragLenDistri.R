# # # # # # # # # # # # # # # # # # # # # # # # # 
#
# script for plot fragment length distribution
#
# author: Weizhang
#
# # # # # # # # # # # # # # # # # # # # # # # # # 

library(rtracklayer)

Args <- commandArgs()
fragfile <- Args[6]
figurefile <- Args[7]

gr <- import(con = fragfile, format = "BED")

frag_len <- end(gr) - start(gr)
frag_freq <- as.data.frame(table(frag_len))
colnames(frag_freq) <- c("Length", "Frequency")

bmp(figurefile)
plot(x = as.numeric(frag_freq$Length), 
     y = as.numeric(frag_freq$Frequency),
     type = "l", lty = 1,
     xlim = c(0 ,1000),
     main = "Fragment Length Distribution",
     xlab = "Fragment Length (bp)",
     ylab = "Frequency")
dev.off()

