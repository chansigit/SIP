# -*- coding: utf-8 -*-
"""
@Time    : 2018/4/3 14:46
@Author  : Weizhang
@FileName: FlowATACAnalysis.py
"""

from ..core import Flow, Report, Configure, Schedule
from ..steps import VarAndClustering


class FlowATACAnalysis(Flow):
    def __init__(self,
                 bamInput=None,
                 peakInput=None,
                 refdir=None,  # for bowtie2, useless in this flow
                 genome=None,  # for motif analysis in chromVAR
                 resultDir='./result',
                 threads=1,
                 ):
        super(FlowATACAnalysis, self).__init__(resultDir=resultDir,
                                               refdir=refdir,
                                               genome=genome,
                                               threads=threads)
        self._setParam('bamInput', bamInput)
        self._setParam('peakInput', peakInput)


    def _call(self,*args):
        bamInput = args[0]._getObj('CellExtracterBam').getOutput('bamOutput')
        peakInput = args[0]._getObj('GenPeakWithFilter').getOutput('bedOutput')
        peakInput = peakInput[0]

        self._setParam('bamInput', bamInput)
        self._setParam('peakInput', peakInput)

    def _build(self,):
        bamInput = self._getParam('bamInput')
        peakInput = self._getParam('peakInput')

        result = VarAndClustering(bamInput=bamInput,
                                  peakInput=peakInput,
                                  genome=self.genome,
                                  threads=self.threads)

        rp = Report()
        rp.add('Section for Variance and Clustering Analysis', [result])

        self._setObj('VarAndClustering', result)


    def _copy(self, ):
        self._linkRecursive(self._getObj('VarAndClustering').getOutput('var.tiff'),
                            self.getFinalRsPath('VarAndClustering'))
        self._linkRecursive(self._getObj('VarAndClustering').getOutput('clustring.tiff'),
                            self.getFinalRsPath('VarAndClustering'))
        self._linkRecursive(self._getObj('VarAndClustering').getOutput('dev.matrix'),
                            self.getFinalRsPath('VarAndClustering'))
        self._linkRecursive(self._getObj('VarAndClustering').getOutput('var.matrix'),
                            self.getFinalRsPath('VarAndClustering'))
        self._linkRecursive(self._getObj('VarAndClustering').getOutput('clu.matrix'),
                            self.getFinalRsPath('VarAndClustering'))








