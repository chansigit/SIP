# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/26 20:54
@Author  : Weizhang
@FileName: FragInPeak.py

count percentage of fragments in the peak region, this
program is a part of quality control for scATAC-seq
"""

from ..core import Step, Configure
import os.path


class FragInPeak(Step):
    def __init__(self,
                 fragInput=None,
                 peakInput=None,
                 reportOutputDir=None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('fragInput', fragInput)
        self.setParamIO('peakInput', peakInput)
        self.setParamIO('reportOutputDir', reportOutputDir)
        self.initIO()

        self._setUpstreamSize(2)

    def impInitIO(self):
        fragInput = self.getParamIO('fragInput')
        peakInput = self.getParamIO('peakInput')
        reportOutputDir = self.getParamIO('reportOutputDir')

        # FragInPeak.R
        self.setInputRscript('FragInPeakR', 'FragInPeak.R')

        # set all input files
        self.setInputDirOrFile('fragInput', fragInput)
        self.setInputDirOrFile('peakInput', peakInput)
        # set output files
        self.setOutputDirNTo1('percentage', None, 'FragInPeak.txt', 'fragInput')

        if reportOutputDir is None:
            self.setParamIO('reportOutputDir', Configure.getTmpDir())

        if peakInput is not None:
            self._setInputSize(len(self.getInputList('peakInput')))

    def call(self, *args):
        fragInput = args[0]
        peakInput = args[1]

        self.setParamIO('fragInput', fragInput.getOutputList('bedOutput'))
        self.setParamIO('peakInput', peakInput.getOutputList('bedOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        fragInput = self.getInputList('fragInput')
        peakInput = self.getInputList('peakInput')
        rScript = self.getInput('FragInPeakR')
        percentage = self.getOutputList('percentage')

        cmdline = [
            'Rscript', rScript,
            os.path.dirname(fragInput[i]),
            peakInput[i], percentage[i]
        ]
        result = self.callCmdline('V1', cmdline)

    def getMarkdownEN(self, ):
        mdtext = """"""
        return mdtext
