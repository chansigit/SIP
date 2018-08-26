# -*- coding: utf-8 -*-
"""
Created on 7th Mar

@author: CyLiu
"""

from ..core import Step, Configure
import subprocess
import os

class TagBarcode(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutput = None,
                 sumOutput = None,
                 baseStart = 1,
                 baseEnd = 12,
                 baseQuality = 10,
                 barcodedRead = 1,
                 discardRead = False,
                 tagName = None,
                 numBaseBelowQuality = 1,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutput', bamOutput)
        self.setParamIO('sumOutput', sumOutput)

        self.setParam('baseStart', baseStart)
        self.setParam('baseEnd', baseEnd)
        self.setParam('baseQuality', baseQuality)
        self.setParam('barcodedRead', barcodedRead)
        self.setParam('discardRead', discardRead)
        self.setParam('tagName', tagName)
        if tagName == 'XC':
            self.setParam('tagType', 'Cell')
        else:
            self.setParam('tagType', 'UMI')
        self.setParam('numBaseBelowQuality', numBaseBelowQuality)
        self.initIO()

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutput = self.getParamIO('bamOutput')
        sumOutput = self.getParamIO('sumOutput')
        tagName = self.getParam('tagName')
        #print(tagName)
        self.setInputDirOrFile('bamInput', bamInput)
        self.setOutputDirNTo1('bamOutput', bamOutput, 'unalign_tagged_%s.bam'%(tagName), 'bamInput')
        self.setOutputDirNTo1('sumOutput', sumOutput, 'unalign_tagged_%s.bam_summary.txt'%(tagName), 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath('unalign_tagged_%s.bam'%(tagName)))
        if sumOutput is None:
            self.setParamIO('sumOutput', Configure.getTmpPath('unalign_tagged_%s.bam_summary.txt'%(tagName)))


    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        bamOutput = self.getOutputList('bamOutput')
        sumOutput = self.getOutputList('sumOutput')

        baseStart = self.getParam('baseStart')
        baseEnd = self.getParam('baseEnd')
        baseQuality = self.getParam('baseQuality')
        barcodedRead = self.getParam('barcodedRead')
        discardRead = self.getParam('discardRead')
        tagName = self.getParam('tagName')
        numBaseBelowQuality = self.getParam('numBaseBelowQuality')

        cmdline = [
                #'/root/software/Drop-seq_tools-1.13/
                'TagBamWithReadSequenceExtended',
                'INPUT=%s'%(bamInput[i]), 'OUTPUT=%s'%(bamOutput[i]), 'SUMMARY=%s'%(sumOutput[i]),
                'BASE_RANGE=%d-%d'%(baseStart, baseEnd), 'BASE_QUALITY=%d'%(baseQuality),
                'BARCODED_READ=%d'%(barcodedRead), 'DISCARD_READ=%s'%(str(discardRead)),
                'TAG_NAME=%s'%(tagName), 'NUM_BASES_BELOW_QUALITY=%d'%(numBaseBelowQuality)
                ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## Tag {tagType} Barcode Result
The {tagType} barcode's quality report is as follow:
```{{r echo=FALSE}}
read.table('{sumOutput}', header=TRUE)
```

        """.format(tagType=self.getParam("tagType"), sumOutput=self.getOutput('sumOutput'))
        return mdtext
