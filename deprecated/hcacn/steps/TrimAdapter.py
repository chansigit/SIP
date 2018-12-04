# -*- coding: utf-8 -*-
"""
Created on 8th Mar

@author: CyLiu
"""

from ..core import Step, Configure
import subprocess
import os

class TrimAdapter(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutput = None,
                 sumOutput = None,
                 adapterSeq = 'AAGCAGTGGTATCAACGCAGAGTACATGGG',
                 #adapterSeq = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
                 misMatches = 0,
                 numBases = 5,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutput', bamOutput)
        self.setParamIO('sumOutput', sumOutput)

        self.setParam('adapterSeq', adapterSeq)
        self.setParam('misMatches', misMatches)
        self.setParam('numBases', numBases)

        self.initIO()

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutput = self.getParamIO('bamOutput')
        sumOutput = self.getParamIO('sumOutput')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('bamOutput', bamOutput, 'unaligned_tagged_trimmed_smart.bam', 'bamInput')
        self.setOutputDirNTo1('sumOutput', sumOutput, 'adapter_trimming_report.txt', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))
        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath('unaligned_tagged_trimmed_smart.bam'))
        if sumOutput is None:
            self.setParamIO('sumOutput', Configure.getTmpPath('adapter_trimming_report.txt'))


    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')
        sumOutput = self.getOutputList('sumOutput')

        adapterSeq = self.getParam('adapterSeq')
        misMatches = self.getParam('misMatches')
        numBases = self.getParam('numBases')

        cmdline = [
                'TrimStartingSequence',
                'INPUT=%s'%(bamInput[i]), 'OUTPUT=%s'%(bamOutput[i]),
                'OUTPUT_SUMMARY=%s'%(sumOutput[i]), 'SEQUENCE=%s'%(adapterSeq),
                'MISMATCHES=%d'%(misMatches), 'NUM_BASES=%d'%(numBases)
        ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## TrimAdapter Result
TrimAdapter's summary report is as follow:
```{{r echo=FALSE}}
library(ggplot2)
taSum <-  file('{sumOutput}','r',blocking=FALSE)
lines <- readLines(taSum)
rowNum <- length(lines)
i <- 1
while(lines[i]=="" || strsplit(lines[i], split="\\t")[[1]][1] != "BIN"){{ i <- i + 1 }}
strlist <- strsplit(lines[(i+1):(rowNum-1)], split="\\t")
strmtx <- do.call(rbind, strlist)
data <- data.frame(BIN=as.integer(strmtx[,1]), VALUE=as.integer(strmtx[,2]))
ggplot(data, aes(x=BIN, y=VALUE)) +
geom_bar(stat = "identity")+
labs( x="Matched Adapter Bins",y="Reads Count" )+theme_bw()+
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
