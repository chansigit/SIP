# -*- coding: utf-8 -*-
"""

@author: Weizhang

merge reads file to fragment file, for merged sample and every single sample

"""

from ..core import Step
from .BamToBed import BamToBed
from .RmChrOrMergeAllSample import RmChrOrMergeAllSample
import re


class MergeToFrag(Step):
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
        self.setOutputDir1To1('bedOutput', bedOutputDir, None, 'frag.bed', 'bedInput')  # do not change!!!

        if bedInput is not None:
            self._setInputSize(len(self.getInputList('bedInput')))

    def call(self, *args):
        samUpstream = args[0]

        if isinstance(samUpstream, BamToBed):
            self.setParamIO('bedInput', samUpstream.getOutput('bedOutput'))
        elif isinstance(samUpstream, RmChrOrMergeAllSample):
            mergedfile = samUpstream.getOutputList('mergedfilename')
            # print("RmChrOrMergeAllSample mergedfilename:")
            # print(mergedfile)
            otherbed = samUpstream.getOutputList('bedOutput')
            # print("RmChrOrMergeAllSample bedOutput:")
            # print(otherbed)
            otherbed = mergedfile + otherbed
            self.setParamIO('bedInput', otherbed)


    def _multiRun(self, ):
        pass

    def _singleRun(self, i):
        bedInput = self.convertToRealPath(self.getInputList('bedInput'))
        bedOutput = self.convertToRealPath(self.getOutputList('bedOutput'))

        tmpfile1 = "file.tmp1.tmp"
        tmpfile2 = "file.tmp2.tmp"
        tmpfile3 = "file.tmp3.tmp"

        f = open(r'%s' % bedInput[i], 'r')
        g = open(r'%s' % tmpfile1, 'w')

        for line in f.readlines():
            temp = re.split(r'\s|/', line)
            new_line = '\t'.join([str(i) for i in temp])
            g.write(new_line + '\n')

        f.close()
        g.close()

        # sort
        cmdline = [
            'sort -k1,1 -k4,4 -k7,7',
            tmpfile1, '>', tmpfile2
        ]
        result = self.callCmdline(None, cmdline)

        # merge
        f = open(r'%s' % tmpfile2, 'r')
        g = open(r'%s' % tmpfile3, 'w')
        for line in f.readlines():
            tmp = line.split()
            if tmp[6] == "-":
                max_v = tmp[2]
            elif tmp[6] == "+":
                min_v = tmp[1]
                if max_v > min_v:
                    new_line = tmp[0] + "\t" + min_v + "\t" + max_v + "\t" + tmp[3] + "\t255\t*\n"
                    g.write(new_line)

        f.close()
        g.close()

        # sort
        cmdline = [
            'sort -k1,1 -k2n,2',
            tmpfile3, '>', bedOutput[i]
        ]
        result = self.callCmdline(None, cmdline)

        # delete tmp file
        cmdline = [
            'rm -f', tmpfile1
        ]
        result = self.callCmdline(None, cmdline)
        cmdline = [
            'rm -f', tmpfile2
        ]
        result = self.callCmdline(None, cmdline)
        cmdline = [
            'rm -f', tmpfile3
        ]
        result = self.callCmdline(None, cmdline)

    def getMarkdownEN(self, ):
        mdtext = """"""
        return mdtext
