# -*- coding: utf-8 -*-
"""
Created on 8th Mar

@author: CyLiu
"""

from ..core import Step, Configure
import subprocess
import os

class BamToFastq(Step):
    def __init__(self,
                 bamInput = None,
                 fastqOutput = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('fastqOutput', fastqOutput)

        self.initIO()

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        fastqOutput = self.getParamIO('fastqOutput')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('fastqOutput', fastqOutput, 'unaligned_mc_tagged_polyA_filtered.fastq', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))
        if fastqOutput is None:
            self.setParamIO('fastqOutput', Configure.getTmpPath('unaligned_mc_tagged_polyA_filtered.fastq'))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        fastqOutput = self.getOutputList('fastqOutput')

        cmdline = [
                'picardBTF %s %s %s'%('Xmx4g', bamInput[i], fastqOutput[i])
        ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
        ## BamToFastq Result
        BamToFastq
        """
        return ""
