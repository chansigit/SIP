# -*- coding: utf-8 -*-
"""
Created on 8th Mar

@author: CyLiu
"""

from ..core import Step, Configure
import subprocess
import os

class FilterBam(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutput = None,
                 tagReject = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutput', bamOutput)

        self.initIO()

        self.setParam('tagReject', tagReject)

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutput = self.getParamIO('bamOutput')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('bamOutput', bamOutput, 'unalign_tagged_filterd.bam', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath('unalign_tagged_filterd.bam'))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')

        tagReject = self.getParam('tagReject')

        cmdline = [
                'FilterBAM', 'INPUT=%s'%(bamInput[i]),
                'OUTPUT=%s'%(bamOutput[i]), 'TAG_REJECT=%s'%(tagReject)
                ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
        ## FilterBam Result
        Discard reads of low quality.
        """
        return ""
