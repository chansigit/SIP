# -*- coding: utf-8 -*-
"""
Created on Six MAR

@author: CyLiu
"""

from ..core import Step,Configure
import subprocess
import os

class FastqToBam(Step):
    def __init__(self,
                 fastqInput1 = None,
                 fastqInput2 = None,
                 bamOutputDir = None,
                 #mapRsOutputDir = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('fastqInput1', fastqInput1)
        self.setParamIO('fastqInput2', fastqInput2)
        self.setParamIO('bamOutputDir', bamOutputDir)
        #self.setParamIO('mapRsOutputDir', mapRsOutputDir)

        self.initIO()


    def impInitIO(self,):
        fastqInput1 = self.getParamIO('fastqInput1')
        fastqInput2 = self.getParamIO('fastqInput2')
        bamOutputDir = self.getParamIO('bamOutputDir')
        #mapRsOutputDir = self.getParamIO('mapRsOutputDir')

        self.setInputDirOrFile('fastqInput1', fastqInput1)
        self.setInputDirOrFile('fastqInput2', fastqInput2)

        self.setOutputDir1To1('bamOutput', bamOutputDir, 'sample', 'bam', 'fastqInput1')
        #self.setOutputDir1To1('mapRsOutput', mapRsOutputDir, 'sample', 'txt', 'fastqInput1')
        if fastqInput1 is not None:
            self._setInputSize(len(self.getInputList('fastqInput1')))

        if bamOutputDir is None:
            self.setParamIO('bamOutputDir', Configure.getTmpDir())

    def call(self, *args):
        fastqUpstream = args[0]
        self.setParamIO('fastqInput1', fastqUpstream.getOutput('fastqOutput1'))
        self.setParamIO('fastqInput2', fastqUpstream.getOutput('fastqOutput2'))

    def _singleRun(self, i):
        fastqInput1 = self.getInputList('fastqInput1')
        fastqInput2 = self.getInputList('fastqInput2')
        bamOutput = self.getOutputList('bamOutput')
        #mapRsOutput = self.getOutputList('mapRsOutput')
        cmdline = ['picardFTB %s %s %s %s'%('Xmx4g',fastqInput1[i], fastqInput2[i], bamOutput[i])]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext = """
## FastqToBam Result
Totally input {num} samples (paired fastq files).
\r\n
        """.format(num = len(self.getInputList('fastqInput1')))
        return mdtext
