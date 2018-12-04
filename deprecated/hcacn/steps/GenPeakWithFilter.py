# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 10:58
@Author  : Weizhang
@FileName: GenPeakWithFilter.py

generate peak from macs2 summit file
summitInput must be a summit file

if topPeak > the number of peaks, R will cause a
"subscript contains out-of-bounds indices" error,
docker report "Execution halted"
"""

from ..core import Step, Configure


class GenPeakWithFilter(Step):
    def __init__(self,
                 summitInput=None,
                 blacklist=None,
                 bedOutputDir=None,
                 overlapRate=0.2,
                 extendRange=250,
                 topPeak=50000,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('summitInput', summitInput)
        self.setParamIO('bedOutputDir', bedOutputDir)
        self.setParamIO('blacklist', blacklist)

        self.initIO()

        # set other parameters
        self.setParam('overlapRate', overlapRate)
        self.setParam('extendRange', extendRange)
        self.setParam('topPeak', topPeak)

    def impInitIO(self):
        summitInput = self.getParamIO('summitInput')
        blacklist = self.getParamIO('blacklist')
        bedOutputDir = self.getParamIO('bedOutputDir')

        # PeakFilter.R
        self.setInputRscript('PeakFilterR', 'PeakFilter.R')

        # set all input files
        self.setInputDirOrFile('summitInput', summitInput)
        self.setInputDirOrFile('blacklist', blacklist)
        # set all output files
        self.setOutputDir1To1('bedOutput', bedOutputDir, None, '_filterd.bed', 'summitInput', '')

        if summitInput is not None:
            self._setInputSize(len(self.getInputList('summitInput')))

    def call(self, *args):
        summitUpstream = args[0]

        self.setParamIO('summitInput', summitUpstream.getOutput('outputSummit'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        summitInput = self.getInputList('summitInput')
        blacklist = self.getInputList('blacklist')
        rScript = self.getInput('PeakFilterR')
        bedOutput = self.getOutputList('bedOutput')

        cmdline = [
            'Rscript', rScript, summitInput[i],
            blacklist[i], bedOutput[i],
            str(self.getParam('overlapRate')), str(self.getParam('extendRange')),
            str(self.getParam('topPeak'))
        ]

        result = self.callCmdline('V1', cmdline)


    def getMarkdownEN(self,):
        summitInput = self.getInputList('summitInput')
        bedOutput = self.getOutputList('bedOutput')
        mdtext ="""
## Peak Filter Result
```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
library(rtracklayer)
summit <- import(con = "{summitInput}", format = "BED")
bed <- import(con = "{bedOutput}", format = "BED")
```
The input peak file contains `r length(summit)` peaks, after filtering, top `r length(bed)` peaks are selected.

        """.format(summitInput=summitInput[0], bedOutput=bedOutput[0])
        return mdtext
