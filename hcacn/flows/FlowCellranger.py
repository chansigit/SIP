# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 11:12:29 2018

@author: Frankie
@date: 20180328
"""
from ..core import Flow, Report, Configure, Schedule
from ..steps import Cellranger, Seuratpreprocessing, Seuratrun
import os

class FlowCellranger(Flow):
    def __init__(self,
                 fastqInput,
                 refdir,
                 genome,
                 resultDir='./result',
                 expectcells=None):
        super(FlowCellranger, self).__init__(resultDir=resultDir,
                                          refdir=refdir,
                                          genome=genome)
        self._setParam('fastqInput', fastqInput)
        self._setParam('expectcells', expectcells)

    def _call(self, *args):
        # args[0]._getObj('FastqDump').getOutput('fastqOutput')
        pass

    def _build(self, ):
        fastqInput = self._getParam('fastqInput')
        expectcells = self._getParam('expectcells')

        test1 = Cellranger(fastqInput = fastqInput, expectcells=expectcells)
        test2 = Seuratpreprocessing()(test1)
        test3 = Seuratrun()(test2)

        rp = Report()
        rp.add('Section for Cellranger', [test1])
        rp.add('Section for Seuratpreprocessing', [test2])
        rp.add('Section for Seuratrun', [test3])

        self._setObj('Cellranger', test1)
        self._setObj('Seuratpreprocessing', test2)
        self._setObj('Seuratrun', test3)
        self._setObj('Report', rp)

    def _copy(self, ):
        self._linkRecursive(self._getObj('Cellranger').getOutput('summary'),self.getFinalRsPath('rawdata_summary'))
        # self._linkRecursive(self._getObj('Cellranger').getOutput('summary'), self.getFinalRsPath('rawdata_summary'))
        # self._linkRecursive(self._getObj('Cellranger').getOutput('summary'), self.getFinalRsPath('rawdata_summary'))
        # self._linkRecursive(self._getObj('Cellranger').getOutput('summary'), self.getFinalRsPath('rawdata_summary'))
        # self._linkRecursive(self._getObj('Cellranger').getOutput('summary'), self.getFinalRsPath('rawdata_summary'))

