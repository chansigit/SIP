rm(list = ls())
library(scde)
args<-commandArgs(T)
gene_expression <- read.table(file = args[1])

cd <- as.matrix(gene_expression)
cd <- apply(cd,2,as.integer)
cd <- clean.counts(cd, min.lib.size=1000, min.reads = 1, min.detected = 1)

sg <- factor(gsub("(cuffquant1|cuffquant2).*", "\\1", colnames(cd)), levels = c("cuffquant1", "cuffquant2"))
names(sg) <- colnames(cd)  
table(sg)
# EVALUATION NOT NEEDED FOR SAKE OF TIME
o.ifm <- scde.error.models(counts = cd, groups = sg, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)

valid.cells <- o.ifm$corr.a > 0

o.ifm <- o.ifm[valid.cells, ]

o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)

groups <- factor(gsub("(cuffquant1|cuffquant2).*", "\\1", rownames(o.ifm)), levels  =  c("cuffquant1", "cuffquant2"))
names(groups) <- row.names(o.ifm)

ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = args[2] + "results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
