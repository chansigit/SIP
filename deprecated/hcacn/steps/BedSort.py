# -*- coding: utf-8 -*-
"""

designed for scATAC-seq

@author: Weizhang

sort BED file using bedtools

"""

from ..core import Step, Configure


class BedSort(Step):
    def __init__(self,
                 bedInput=None,
                 bedOutputDir=None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bedInput', bedInput)
        self.setParamIO('bedOutputDir', bedOutputDir)

        self.initIO()

    def impInitIO(self):
        bedInput = self.getParamIO('bedInput')
        bedOutputDir = self.getParamIO('bedOutputDir')

        # set all input files
        self.setInputDirOrFile('bedInput', bedInput)
        # set all output files
        self.setOutputDir1To1('bedOutput', bedOutputDir, None, 'bed', 'bedInput')

        if bedInput is not None:
            self._setInputSize(len(self.getInputList('bedInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('bedInput', samUpstream.getOutput('mergedfilename'))

    def _multiRun(self,):
        pass
        # bedInput = self.getInputList('bedInput')
        # bedOutput = self.getOutputList('bedOutput')
        # for i in range(len(bedInput)):
        #     cmdline = [
        #         'sortBed -i',
        #         bedInput[i],
        #         '>',
        #         bedOutput[i]
        #     ]
        #     result = self.callCmdline(cmdline)

    def _singleRun(self, i):
        bedInput = self.getInputList('bedInput')
        bedOutput = self.getOutputList('bedOutput')
        cmdline = [
            'sortBed -i',
            bedInput[i],
            '>',
            bedOutput[i]
        ]
        result = self.callCmdline('V1', cmdline)

    def getMarkdownEN(self, ):
        mdtext = """"""
        return mdtext
