# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/29 8:39
@Author  : Weizhang
@FileName: CellFilter.py
"""


from ..core import Step, Configure
import os.path

class CellFilter(Step):
    def __init__(self,
                 libCpxInput=None,
                 fragInPeakInput=None,
                 filterOutputDir=None,
                 libSizeCutOff=10000,  # number of reads from picardLibComplexity
                 fragInPeakCutOff=0.15,  # 15% reads in peak region
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('libCpxInput', libCpxInput)
        self.setParamIO('fragInPeakInput', fragInPeakInput)
        self.setParamIO('filterOutputDir', filterOutputDir)
        self.initIO()

        # set other parameters
        self.setParam('libSizeCutOff', libSizeCutOff)
        self.setParam('fragInPeakCutOff', fragInPeakCutOff)

        self._setUpstreamSize(2)

    def impInitIO(self):
        libCpxInput = self.getParamIO('libCpxInput')
        fragInPeakInput = self.getParamIO('fragInPeakInput')
        filterOutputDir = self.getParamIO('filterOutputDir')

        # CellFilter.R
        self.setInputRscript('CellFilterR', 'CellFilter.R')

        # set all input files
        self.setInputDirOrFile('libCpxInput', libCpxInput)
        self.setInputDirOrFile('fragInPeakInput', fragInPeakInput)
        # set all output files
        self.setOutputDirNTo1('reportOutput', None, 'CellFilter.txt', 'libCpxInput')
        self.setOutputDirNTo1('figureOutput', None, 'CellFilter.bmp', 'libCpxInput')

        if filterOutputDir is None:
            self.setParamIO('filterOutputDir', Configure.getTmpDir())

        if fragInPeakInput is not None:
            self._setInputSize(len(self.getInputList('fragInPeakInput')))

    def call(self, *args):
        libCpxUpstream = args[0]
        fragInPeakUpstream = args[1]

        # samOutput is from the former step (Mapping)
        self.setParamIO('libCpxInput', libCpxUpstream.getOutput('sumOutput'))
        self.setParamIO('fragInPeakInput', fragInPeakUpstream.getOutput('percentage'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        libCpxInput = self.getInputList('libCpxInput')  # list
        fragInPeakInput = self.getInputList('fragInPeakInput')  # only one file
        reportOutput = self.getOutputList('reportOutput')  # a file name
        figureOutput = self.getOutputList('figureOutput')  # a file name
        rScript = self.getInput('CellFilterR')

        cmdline = [
            'Rscript', rScript,
            os.path.dirname(libCpxInput[i]), fragInPeakInput[i],
            reportOutput[i], figureOutput[i],
            str(self.getParam('libSizeCutOff')), str(self.getParam('fragInPeakCutOff'))
        ]
        result = self.callCmdline('V1', cmdline)

    def getMarkdownEN(self, ):

        figureOutput = self.getOutputList('figureOutput')

        mdtext = """
## Cell Filtering Result

The following is a figure about library complexity and reads in peak regions.

The cells located in top right region are selected for downstream analysis.
```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
figure_path <- "{figureOutput}"
```

![](`r figure_path`)

        """.format(figureOutput=figureOutput[0])
        return mdtext
