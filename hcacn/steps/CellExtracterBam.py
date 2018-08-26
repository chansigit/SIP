# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/29 15:12
@Author  : Weizhang
@FileName: CellExtracterBam.py

copy files

"""

from ..core import Step, Configure
import os.path

# x is a str/substr, y is a list
def strin(x, y):
    for tmp_tsr in y:
        if x in tmp_tsr:
            return tmp_tsr

    return 0

class NoFilterdBAM(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class CellExtracterBam(Step):
    def __init__(self,
                 filterFile=None,
                 bamInput=None,
                 outputDir=None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('filterFile', filterFile)
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('outputDir', outputDir)
        self.initIO()

        self._setUpstreamSize(2)

    def impInitIO(self):
        filterFile = self.getParamIO('filterFile')
        bamInput = self.getParamIO('bamInput')
        outputDir = self.getParamIO('outputDir')

        # set all input files
        self.setInputDirOrFile('filterFile', filterFile)
        self.setInputDirOrFile('bamInput', bamInput)
        # set all output files
        self.setOutputDir1To1('bamOutput', outputDir, None, 'bam', 'bamInput')

        if outputDir is None:
            self.setParamIO('outputDir', Configure.getTmpDir())

        if filterFile is not None:  # flag1
            self._setInputSize(len(self.getInputList('filterFile')))

    def _preRun(self,):
        filter_file = self.getInputList('filterFile')  # a str, file name
        filter_file = filter_file[0]
        bam_file = self.getInputList('bamInput')  # a list, list of file name
        save_file = []
        outputDir = self.getParamIO('outputDir')

        flag_for_empty = os.path.getsize(filter_file)
        if flag_for_empty == 0:
            raise NoFilterdBAM('No file could pass the CellFilter step, please see the figure from '
                               'CellFilter step and change your threshold!!')

        f = open(filter_file, 'r')
        for line in f.readlines():
            filenameflag = line.split()[0]
            tmp_saved = strin(filenameflag, bam_file)
            save_file.append(tmp_saved)

        f.close()
        self.setInputDirOrFile('bamInput', save_file)

        # set all output files
        self.setOutputDir1To1('bamOutput', outputDir, None, 'bam', 'bamInput')

        if outputDir is None:
            self.setParamIO('outputDir', Configure.getTmpDir())

        if save_file is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

    def call(self, *args):
        filterUpstream = args[0]
        bamUpstream = args[1]

        # samOutput is from the former step (Mapping)
        self.setParamIO('filterFile', filterUpstream.getOutput('reportOutput'))
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')

        cmdline = [
            'ln', '-f', bamInput[i], bamOutput[i]
        ]
        result = self.callCmdline('V1', cmdline, stdoutToLog=False)

    def getMarkdownEN(self, ):
        mdtext = """"""
        return mdtext
