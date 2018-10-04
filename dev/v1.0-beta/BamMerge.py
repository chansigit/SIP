

from ..core import Step,Configure
import subprocess
import os

class BamMerge(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutput = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutput', bamOutput)

        self.initIO()

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutput = self.getParamIO('bamOutput')

        self.setInputDirOrFile('bamInput', bamInput)
        #self.setOutputDir1To1('bamOutput', bamOutputDir, 'merged', 'bam', 'bamInput')
        self.setOutputDirNTo1('bamOutput', bamOutput, 'unaligned_data.bam', 'bamInput')
        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath('unaligned_data.bam'))

        self._setMultiRun()

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _multiRun(self,):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')
        cmdline = ['samtools merge', bamOutput[0]]
        for i in range(len(bamInput)):
            cmdline.append(bamInput[i])
        result = self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
        ## BamMerge Result
        BamMerge
        """
        return ""
