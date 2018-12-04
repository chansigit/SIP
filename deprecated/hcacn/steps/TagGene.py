# -*- coding: utf-8 -*-
"""
Created on 9th Mar

@author: CyLiu
"""

from ..core import Step, Configure
import subprocess
import os

class TagGene(Step):
    def __init__(self,
                 bamInput = None,
                 bamOutput = None,
                 gtfInput = None,
                 tag = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        Configure.setRefSuffix('gtfRef', '_refdata-cellranger-2.1.0/genes/genes.gtf', check=False)
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('bamOutput', bamOutput)
        self.setParamIO('gtfInput', gtfInput)

        self.setParam('tag', tag)

        self.initIO()

    def impInitIO(self,):
        bamInput = self.getParamIO('bamInput')
        bamOutput = self.getParamIO('bamOutput')
        gtfInput = self.getParamIO('gtfInput')
        if gtfInput is None:
            self.setParamIO('gtfInput', Configure.getConfig('gtfRef'))
            gtfInput = self.getParamIO('gtfInput')

        self.setInputDirOrFile('bamInput', bamInput)
        self.setInputDirOrFile('gtfInput', gtfInput)
        self.setOutputDirNTo1('bamOutput', bamOutput, 'star_gene_exon_tagged.bam', 'bamInput')

        if bamInput is not None:
            self._setInputSize(len(self.getInputList('bamInput')))

        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath('star_gene_exon_tagged.bam'))

    def call(self, *args):
        bamUpstream = args[0]
        self.setParamIO('bamInput', bamUpstream.getOutput('bamOutput'))

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        gtfInput = self.getInputList('gtfInput')
        bamOutput = self.getOutputList('bamOutput')

        tag = self.getParam('tag')

        cmdline = [
                'TagReadWithGeneExon', 'I=%s'%(bamInput[i]),
                'O=%s'%(bamOutput[i]), 'ANNOTATIONS_FILE=%s'%(gtfInput[i]), 'TAG=%s'%(tag)
        ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
        ## TagGene Result
        Tag reads with gene exon.
        """
        return ""
