# -*- coding: utf-8 -*-
from ..core import Step, Configure
import subprocess
import os

class DetectError(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutput = None,
                 statsOutput = None,
                 sumOutput = None,
                 numCells = 1000,
                 primerSeqence = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutput', bamOutput)
        self.setParamIO('statsOutput', statsOutput)
        self.setParamIO('sumOutput', sumOutput)

        self.initIO()

        self.setParam('numCells', numCells)
        self.setParam('primerSeqence', primerSeqence)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutput = self.getParamIO('bamOutput')
        statsOutput = self.getParamIO('statsOutput')
        sumOutput = self.getParamIO('sumOutput')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('bamOutput', bamOutput, 'out_gene_exon_tagged.bam', 'bamInput')
        self.setOutputDirNTo1('statsOutput', statsOutput, 'synthesis_stats.txt', 'bamInput')
        self.setOutputDirNTo1('sumOutput', sumOutput, 'synthesis_stats.summary.txt', 'bamInput')
        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))
        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath('out_gene_exon_tagged.bam'))
        if statsOutput is None:
            self.setParamIO('statsOutput', Configure.getTmpPath('synthesis_stats.txt'))
        if sumOutput is None:
            self.setParamIO('sumOutput', Configure.getTmpPath('synthesis_stats.summary.txt'))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')
        statsOutput = self.getOutputList('statsOutput')
        sumOutput = self.getOutputList('sumOutput')

        numCells = self.getParam('numCells')
        primerSeqence = self.getParam('primerSeqence')

        cmdline = [
                'DetectBeadSynthesisErrors',
                'I=%s'%(bamInput[i]), 'O=%s'%(bamOutput[i]), 'OUTPUT_STATS=%s'%(statsOutput[i]),
                'SUMMARY=%s'%(sumOutput[i]), 'NUM_BARCODES=%d'%(4*numCells),
                'PRIMER_SEQUENCE=%s'%(primerSeqence)
        ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## DetectError Result
Synthesis status summary report is shown below:
```{{r echo=FALSE}}
deSum <- file('{sumOutput}', 'r', blocking=FALSE)
lines <- readLines(deSum)
i <- 1
while(lines[i]=="" || strsplit(lines[i], split="\t")[[1]][1] != "NUM_BEADS"){{ i <- i + 1 }}
strlist <- strsplit(lines[(i):(i+1)], split="\t")
strmtx <- do.call(rbind, strlist)
data.frame(Info=strmtx[1,], Count=as.integer(strmtx[2,]))
```

        """.format(sumOutput=self.getOutput('sumOutput'))
        return mdtext
