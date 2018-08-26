# -*- coding: utf-8 -*-
"""

designed for scATAC-seq

@author: Weizhang

"""

from ..core import Step


class BamToBed(Step):
    def __init__(self,
                 bamInput=None,
                 bedOutputDir=None,
                 FSShift=4,  # forward strand shift
                 RSShift=-5,  # reverse strand shift
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bedOutputDir', bedOutputDir)

        self.initIO()

        # set other parameters
        if FSShift >= 0:
            self.setParam('FSShift', '+ ' + str(FSShift))  # save as a string "+ 4"
        else:
            self.setParam('FSShift', '- ' + str(abs(FSShift)))  # save as a string "- 5"
        if RSShift >= 0:
            self.setParam('RSShift', '+ ' + str(RSShift))  # save as a string "+ 4"
        else:
            self.setParam('RSShift', '- ' + str(abs(RSShift)))  # save as a string "- 5"

    def impInitIO(self):
        bamInput = self.getParamIO('bamInput')
        bedOutputDir = self.getParamIO('bedOutputDir')

        # set all input files
        self.setInputDirOrFile('bamInput', bamInput)
        # set all output files
        self.setOutputDir1To1('bedOutput', bedOutputDir, None, 'bed', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        samUpstream = args[0]

        # samOutput is from the former step (Mapping)
        self.setParamIO('bamInput', samUpstream.getOutput('bamOutput'))

    def _multiRun(self,):
        pass
        # bamInput = self.getInputList('bamInput')
        # bedOutput = self.getOutputList('bedOutput')
        #
        # for i in range(len(bamInput)):
        #     cmdline = [
        #         "bedtools bamtobed -i", bamInput[i],
        #         "|",
        #         "awk \'BEGIN {OFS = \"\\t\"} ; {if ($6 == \"+\") print $1, $2",
        #         self.getParam('FSShift'),  # "+ 4"
        #         ", $3",
        #         self.getParam('FSShift'),  # "+ 4"
        #         ", $4, $5, $6; else print $1, $2",
        #         self.getParam('RSShift'),  # "- 5"
        #         ", $3",
        #         self.getParam('RSShift'),  # "- 5"
        #         ", $4, $5, $6}\' - >",
        #         bedOutput[i]
        #     ]
        #     result = self.callCmdline(cmdline)

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bedOutput = self.getOutputList('bedOutput')

        cmdline = [
            "bedtools bamtobed -i", bamInput[i],
            "|",
            "awk \'BEGIN {OFS = \\\"\\\\t\\\"} ; {if (\$6 == \\\"+\\\") print \$1, \$2",
            self.getParam('FSShift'),  # "+ 4"
            ", \$3",
            self.getParam('FSShift'),  # "+ 4"
            ", \$4, \$5, \$6; else print \$1, \$2",
            self.getParam('RSShift'),  # "- 5"
            ", \$3",
            self.getParam('RSShift'),  # "- 5"
            ", \$4, \$5, \$6}\' - >",
            bedOutput[i]
        ]
        result = self.callCmdline('V1', cmdline)


    def getMarkdownEN(self, ):
        mdtext = """"""
        return mdtext
