# -*- coding: utf-8 -*-
"""

@author: Weizhang

this program convert SAM to BAM, including reads filter

"""

from ..core import Step, Configure
import os


class SamToBam(Step):
    # set default function: convert SAM to BAM
    def __init__(self,
                 samInput=None,  # <in.bam>|<in.sam>|<in.cram>
                 bamOutputDir=None,  # -o FILE  output file name [stdout]
                 threads=1,  # -@ INT number of BAM compression threads [0]
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('samInput', samInput)
        self.setParamIO('bamOutputDir',bamOutputDir)

        self.initIO()

        # set other parameters
        self.setParam('threads', threads)

    def impInitIO(self):
        samInput = self.getParamIO('samInput')
        bamOutputDir = self.getParamIO('bamOutputDir')
        if bamOutputDir is None:
            self.setParamIO('bamOutputDir',Configure.getTmpDir())
            
        self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))

        # set all input files
        self.setInputDirOrFile('samInput', samInput)
        # set all output files
        self.setOutputDir1To1('bamOutput', bamOutputDir, None, 'bam', 'samInput')

        if samInput is not None:
            self._setInputSize(len(self.getInputList('samInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('samInput', samUpstream.getOutput('samOutput'))

    def _multiRun(self,):
        pass
        # samInput = self.getInputList('samInput')
        # bamOutput = self.getOutputList('bamOutput')
        # for i in range(len(samInput)):
        #     cmdline = [
        #         'samtools view -b -S',
        #         '-@', str(self.getParam('threads')),
        #         '-o', bamOutput[i], samInput[i]
        #     ]
        #     result = self.callCmdline(cmdline)

    def _singleRun(self, i):
        samInput = self.getInputList('samInput')
        bamOutput = self.getOutputList('bamOutput')

        cmdline = [
            'samtools view -bS -q 30 -F 1804',
            '-@', str(self.getParam('threads')),
            '-o', bamOutput[i], samInput[i]
        ]

        result = self.callCmdline('V1', cmdline)
        f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
        f.write(result.stdout)
        f.close()
            
    def getMarkdownEN(self,):
        # this function do not need any report!
        mdtext = """"""
        return mdtext