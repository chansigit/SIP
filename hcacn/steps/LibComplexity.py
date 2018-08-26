# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/26 18:02
@Author  : Weizhang
@FileName: LibComplexity.py
"""

from ..core import Step, Configure


class LibComplexity(Step):
    def __init__(self,
                 bamInput=None,
                 sumOutputDir=None,
                 memory='-Xmx5g',
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('sumOutputDir', sumOutputDir)

        self.initIO()

        # set other parameters
        self.setParam('memory', memory)

    def impInitIO(self):
        bamInput = self.getParamIO('bamInput')
        sumOutputDir = self.getParamIO('sumOutputDir')

        # set all input files
        self.setInputDirOrFile('bamInput', bamInput)
        # set all output files
        self.setOutputDir1To1('sumOutput', sumOutputDir, None, 'lcpx', 'bamInput')

        if sumOutputDir is None:
            self.setParamIO('sumOutputDir', Configure.getTmpDir())

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('bamInput', samUpstream.getOutput('bamOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        sumOutput = self.getOutputList('sumOutput')
        cmdline = [
            'picardLibComplexity', str(self.getParam('memory')),
            bamInput[i], sumOutput[i]
        ]
        result = self.callCmdline('V1', cmdline)

    def getMarkdownEN(self, ):
        mdtext = """"""
        return mdtext
