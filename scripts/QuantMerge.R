# Args[6]=={File path list}
# Args[7]=={Output}
# Args[8]=={t2g file}
# Args[9]=={quant unit}
# example:
#     Rscript QuantMerge.R /data/hca/chensijie/SMART1_small_products/Quant.salmon/PathListFile.txt  DGE.tsv /data/hca/chensijie/code/nextflow_test/scripts/tx2gene.SALMON21_Homo_sapiens.GRCh38.cdna.all.tsv  count 

Args <- commandArgs()
filepath <-as.character(Args[6])
output  <-as.character(Args[7])
t2gFile <- as.character(Args[8])
quantunit<-as.character(Args[9])

library(tximport)
library(stringr)

# 1. Load data
##  Load the path list file from command line parameter, and read the path list file
#filepath="/data/hca/chensijie/SMART1_small_products/Quant.salmon/PathListFile.txt"
txt <- readLines(filepath)


# 2. Remove the two square brackets at the beginning/end
## The path list file contains comma separated paths contained in-between left/right square brackets,
## like [/path/to/cell1_quant/quant.sf, /path/to/cell2_quant/quant.sf, /path/to/cell3_quant/quant.sf]
txt     <- gsub("\\[|\\]", "", txt)  # The txt file contains 
fileDF  <- read.csv(text = txt, header = FALSE)
fileVec <- as.vector(t(fileDF))
fileVec <- str_trim(fileVec)
## this produces a character vector like: 
## c('/path/to/cell1_quant/quant.sf', '/path/to/cell2_quant/quant.sf', '/path/to/cell3_quant/quant.sf')


# 3. Set cell names and sort cell by its identifier
# We should now extract cell identifier from paths and set them as vector names.
# Retrieve the parent folder and trim the _quant tail
cell.id   <- function (path) basename(dirname(path))
tail.trim <- function (x) sub("_quant$", "", x)
names(fileVec) <- tail.trim(cell.id(fileVec))
fileVec <- fileVec[order(names(fileVec))]


# 4. Load tx2gene file
## This tx2gene checking table must lie in the same path as this R script
## We may use another tx2gene file when using different annotations.
tx2gene <-   read.table(t2gFile)

## A simple list with matrices,"abundance","counts", and "length", is returned,
## where the transcript level information is summarized to the gene-level.

## The “length” matrix can be used to generate an offset matrix 
## for downstream gene-level differential analysis of count matrices
genelevelSummary <- tximport(fileVec, type="salmon", tx2gene=tx2gene)

## Select the count matrix for output
if (quantunit=="count"){
    DGE    <- genelevelSummary$counts
    write.table(DGE,  file=output,    sep="\t")
}
