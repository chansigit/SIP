# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""
from ..core import Flow,Report,Configure,Schedule
from ..steps import MonocleDC,MonocleQC
      
class MonocleDimentionReductionAndClusterFlow(Flow):
    def __init__(self,
                 matrixdata,
                 min_expression = 0.1,
                 num_cells_expressed_threshold = 10,
                 TotalmRNAs = 1e6, 
                 mean_expression_threshold=0.1,
                 num_PCA = 10,
                 cluster_num = 6,
                 resultDir='./result'):
        super(MonocleDimentionReductionAndClusterFlow,self).__init__(resultDir=resultDir, 
                                                                     refdir=None, 
                                                                     genome=None, 
                                                                     threads=None)
        self._setParam('matrixdata',matrixdata)
        self._setParam('min_expression',min_expression)
        self._setParam('num_cells_expressed_threshold',num_cells_expressed_threshold)
        self._setParam('TotalmRNAs',TotalmRNAs)
        self._setParam('mean_expression_threshold',mean_expression_threshold)
        self._setParam('num_PCA',num_PCA)
        self._setParam('cluster_num',cluster_num)
        
    def _call(self,*args):
        # args[0]._getObj('FastqDump').getOutput('fastqOutput')
        pass
    
    def _build(self,):
        matrixdata = self._getParam('matrixdata')
        min_expression = self._getParam('min_expression')
        num_cells_expressed_threshold = self._getParam('num_cells_expressed_threshold')
        TotalmRNAs = self._getParam('TotalmRNAs')
        mean_expression_threshold = self._getParam('mean_expression_threshold')
        num_PCA = self._getParam('num_PCA')
        cluster_num = self._getParam('cluster_num')
        MonocleQC_result = MonocleQC(matrixdata = matrixdata, 
                                     outputpath = None,
                                     min_expression = min_expression,
                                     num_cells_expressed_threshold = num_cells_expressed_threshold,
                                     TotalmRNAs = TotalmRNAs, 
                                     mean_expression_threshold = mean_expression_threshold)
        MonocleDC_result = MonocleDC(outputpath = None,
                                     num_PCA = num_PCA,
                                     cluster_num = cluster_num)(MonocleQC_result)
        rp = Report()
        rp.add('Monocle Quality Control Results',[MonocleQC_result])
        rp.add('Monocle Dimension Reduction And Density Peak Clustering Results',[MonocleDC_result])
        self._setObj('MonocleQC_result',MonocleQC_result)
        self._setObj('MonocleDC_result',MonocleDC_result)
        
    def _copy(self,):
        pass