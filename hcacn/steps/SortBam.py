# -*- coding: utf-8 -*-
"""
Created on 9th Mar

@author: CyLiu
"""

from ..core import Step, Configure
import subprocess
import os

class SortBam(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutput = None,
                 sortOrder = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutput', bamOutput)
        self.setParam('sortOrder', sortOrder)

        self.initIO()


    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutput = self.getParamIO('bamOutput')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('bamOutput', bamOutput, 'aligned.sorted.bam', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))
        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath('aligned.sorted.bam'))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')

        sortOrder = self.getParam('sortOrder')

        cmdline = [
                'picardSB %s %s %s %s'%('Xmx4g', bamInput[i], bamOutput[i], sortOrder)
        ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
        ## SortBam Result
        Sort bam file as query order.
        """
        return ""
