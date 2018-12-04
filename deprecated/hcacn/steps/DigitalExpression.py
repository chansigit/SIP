# -*- coding: utf-8 -*-

from ..core import Step, Configure
import subprocess
import os

class DigitalExpression(Step):
    def __init__(self,
                 bamInput = None,
                 dgeOutput = None,
                 sumOutput = None,
                 numCells = 1000,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('sumOutput', sumOutput)
        self.setParamIO('dgeOutput', dgeOutput)
        if dgeOutput is not None:
            self.setParamIO('midOutput', dgeOutput+".gz")
        else:
            self.setParamIO('midOutput', None)
        self.setParam('numCells', numCells)

        self.initIO()

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        midOutput = self.getParamIO('midOutput')
        sumOutput = self.getParamIO('sumOutput')
        dgeOutput = self.getParamIO('dgeOutput')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('midOutput', midOutput, 'out_gene_exon_tagged.dge.txt.gz', 'bamInput')
        self.setOutputDirNTo1('sumOutput', sumOutput, 'out_gene_exon_tagged.dge.summary.txt', 'bamInput')
        self.setOutputDirNTo1('dgeOutput', dgeOutput, 'out_gene_exon_tagged.dge.txt', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))
        if midOutput is None:
            self.setParamIO('midOutput', Configure.getTmpPath('out_gene_exon_tagged.dge.txt.gz'))
        if sumOutput is None:
            self.setParamIO('sumOutput', Configure.getTmpPath('out_gene_exon_tagged.dge.summary.txt'))
        if dgeOutput is None:
            self.setParamIO('dgeOutput', Configure.getTmpPath('out_gene_exon_tagged.dge.txt'))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        midOutput = self.getOutputList('midOutput')
        sumOutput = self.getOutputList('sumOutput')
        dgeOutput = self.getOutputList('dgeOutput')

        numCells = self.getParam('numCells')

        cmdline = [
                'DigitalExpression',
                'I=%s'%(bamInput[i]), 'O=%s'%(midOutput[i]), 'SUMMARY=%s'%(sumOutput[i]),
                'NUM_CORE_BARCODES=%d'%(numCells)
        ]
        self.callCmdline('V1', cmdline)

        cmdline = ['gunzip -c %s'%(midOutput[i]), '>', dgeOutput[i]]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## DigitalExpression Result
The summary os digitalExpression is shown below:

### The distribution of number of genic reads per cell:
```{{r echo=FALSE}}
library(ggplot2)
dgeSum <- file('{sumOutput}', "r", blocking=FALSE)
lines <- readLines(dgeSum)
rowNum <- length(lines)
i <- 1
while(lines[i]=="" || strsplit(lines[i], split="\\t")[[1]][1] != "CELL_BARCODE"){{ i <- i + 1 }}
strlist <- strsplit(lines[(i+1):(rowNum-1)], split="\\t")
strmtx <- do.call(rbind, strlist)
#data <- data.frame(CELL_BARCODE=strmtx[,1], NUM_GENIC_READS=as.integer(strmtx[,2]),
#NUM_TRANSCRIPTS=as.integer(strmtx[,3]), NUM_GENES=as.integer(strmtx[,4]))
data.plot <- data.frame(count = as.integer(strmtx[,2]))
ggplot( data.plot,aes(x=count) )+
geom_histogram(bins=100,fill='lightblue1',color='gray',aes(y=..count..))+
labs( x="Number of genic reads",y="Counts" )+theme_bw()+scale_x_log10()+
annotate( "text",label=paste0("Meidan: ",median(data.plot$count)),size=5,colour='black',x=max(data.plot$count)/2,y=15 )+
geom_density( aes(y=0.01*..count..) )+
theme(plot.title=element_text(size=20),
axis.title.y=element_text(size = 16, vjust=+0.2),
axis.title.x=element_text(size = 16, vjust=-0.2),
axis.text.y=element_text(size = 14),
axis.text.x=element_text(size = 14),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank())
```

### The distribution of number of transcripts per cell:
```{{r echo=FALSE}}
library(ggplot2)
dgeSum <- file('{sumOutput}', "r", blocking=FALSE)
lines <- readLines(dgeSum)
rowNum <- length(lines)
i <- 1
while(lines[i]=="" || strsplit(lines[i], split="\\t")[[1]][1] != "CELL_BARCODE"){{ i <- i + 1 }}
strlist <- strsplit(lines[(i+1):(rowNum-1)], split="\\t")
strmtx <- do.call(rbind, strlist)
#data <- data.frame(CELL_BARCODE=strmtx[,1], NUM_GENIC_READS=as.integer(strmtx[,2]),
#NUM_TRANSCRIPTS=as.integer(strmtx[,3]), NUM_GENES=as.integer(strmtx[,4]))
data.plot <- data.frame(count = as.integer(strmtx[,3]))
ggplot( data.plot,aes(x=count) )+
geom_histogram(bins=100,fill='lightblue1',color='gray',aes(y=..count..))+
labs( x="Number of transcripts",y="Counts" )+theme_bw()+scale_x_log10()+
annotate( "text",label=paste0("Meidan: ",median(data.plot$count)),size=5,colour='black',x=max(data.plot$count)/2,y=15 )+
geom_density( aes(y=0.01*..count..) )+
theme(plot.title=element_text(size=20),
axis.title.y=element_text(size = 16, vjust=+0.2),
axis.title.x=element_text(size = 16, vjust=-0.2),
axis.text.y=element_text(size = 14),
axis.text.x=element_text(size = 14),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank())
```

### The distribution of number of genes per cell:
```{{r echo=FALSE}}
library(ggplot2)
dgeSum <- file('{sumOutput}', "r", blocking=FALSE)
lines <- readLines(dgeSum)
rowNum <- length(lines)
i <- 1
while(lines[i]=="" || strsplit(lines[i], split="\\t")[[1]][1] != "CELL_BARCODE"){{ i <- i + 1 }}
strlist <- strsplit(lines[(i+1):(rowNum-1)], split="\\t")
strmtx <- do.call(rbind, strlist)
#data <- data.frame(CELL_BARCODE=strmtx[,1], NUM_GENIC_READS=as.integer(strmtx[,2]),
#NUM_TRANSCRIPTS=as.integer(strmtx[,3]), NUM_GENES=as.integer(strmtx[,4]))
data.plot <- data.frame(count = as.integer(strmtx[,4]))
ggplot( data.plot,aes(x=count) )+
geom_histogram(bins=100,fill='lightblue1',color='gray',aes(y=..count..))+
labs( x="Number of genes",y="Counts" )+theme_bw()+
annotate( "text",label=paste0("Meidan: ",median(data.plot$count)),size=5,colour='black',x=max(data.plot$count)/2,y=15 )+
geom_density( aes(y=30*..count..) )+
theme(plot.title=element_text(size=20),
axis.title.y=element_text(size = 16, vjust=+0.2),
axis.title.x=element_text(size = 16, vjust=-0.2),
axis.text.y=element_text(size = 14),
axis.text.x=element_text(size = 14),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank())
```

        """.format(sumOutput=self.getOutput('sumOutput'))
        return mdtext
