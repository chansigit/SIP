# -*- coding: utf-8 -*-

from ..core import Step, Configure
import subprocess
import os

class MergeBamAlign(Step):
    def __init__(self,
                 unmappedBamInput = None,
                 alignedBamInput = None,
                 bamOutput = None,
                 refInputDir = None,
                 secondAlign = False,
                 pairedRun = False,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)
        Configure.setRefSuffix('mergeRef', '_refdata-cellranger-2.1.0/fasta', check=False)
        self.setParamIO('unmappedBamInput', unmappedBamInput)
        self.setParamIO('alignedBamInput', alignedBamInput)
        self.setParamIO('bamOutput', bamOutput)
        self.setParamIO('refInputDir', refInputDir)
        self.setParam('secondAlign', secondAlign)
        self.setParam('pairedRun', pairedRun)

        self.initIO()


    def impInitIO(self,):
        unmappedBamInput = self.getParamIO('unmappedBamInput')
        alignedBamInput = self.getParamIO('alignedBamInput')
        bamOutput = self.getParamIO('bamOutput')
        refInputDir = self.getParamIO('refInputDir')
        if refInputDir is None:
            self.setParamIO('refInputDir', Configure.getConfig('mergeRef'))
            refInputDir = self.getParamIO('refInputDir')
        self.setInputDirOrFile('unmappedBamInput', unmappedBamInput)
        self.setInputDirOrFile('alignedBamInput', alignedBamInput)
        self.setInputDirOrFile('refSeqInput', os.path.join(refInputDir, 'genome.fa'))
        self.setInputDirOrFile('refDictInput', os.path.join(refInputDir, 'genome.dict'))
        self.setOutputDirNTo1('bamOutput', bamOutput, 'merged.bam', 'unmappedBamInput')

        self._setUpstreamSize(2)

        if unmappedBamInput is not None:
            self._setInputSize(len(self.getInputList('unmappedBamInput')))

        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath('merged.bam'))


    def call(self, *args):
        unmappedBamUpstream = args[0]
        alignedBamUpstream = args[1]

        self.setParamIO('unmappedBamInput', unmappedBamUpstream.getOutput('bamOutput'))
        self.setParamIO('alignedBamInput', alignedBamUpstream.getOutput('bamOutput'))

    def _singleRun(self,i):
        unmappedBamInput = self.getInputList('unmappedBamInput')
        alignedBamInput = self.getInputList('alignedBamInput')
        bamOutput = self.getOutputList('bamOutput')
        refSeqInput = self.getInputList('refSeqInput')

        secondAlign = self.getParam('secondAlign')
        pairedRun = self.getParam('pairedRun')

        cmdline = [
                'picardMBA %s %s %s %s %s %s %s'%('Xmx4g', unmappedBamInput[i], alignedBamInput[i], bamOutput[i], refSeqInput[i],
                str(secondAlign), str(pairedRun))
        ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
        ## MergeBamAlign Result
        Merge unaligned bam and sorted bam.
        """
        return ""
