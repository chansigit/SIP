# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""
from ..core import Flow,Report,Configure,Schedule
from ..steps import MonocleDC,MonocleQC,SingleCellExperiment,SC3_Cluster
      
class Cluster(Flow):
    def __init__(self,
                 #select the algorithm
                 algorithm = 'Monocle',
                 #Monocle parameters list
                 matrixdata = None,
                 min_expression = 0.1,
                 num_cells_expressed_threshold = 10,
                 TotalmRNAs = 1e6, 
                 mean_expression_threshold=0.1,
                 num_PCA = 10,
                 mncluster_num = 6,
                 #sc3 parameters list
                 sc3matrix_file= None,
                 matrix_format = 'ORIGIN',
                 sc3cluster_num = 0,
                 sc3ann = "",
                 resultDir='./result'):
        super(Cluster,self).__init__(resultDir=resultDir, 
                                     refdir=None, 
                                     genome=None, 
                                     threads=None)
        self._setParam('algorithm', algorithm)
        #set Monocle parameters
        self._setParam('matrixdata',matrixdata)
        self._setParam('min_expression',min_expression)
        self._setParam('num_cells_expressed_threshold',num_cells_expressed_threshold)
        self._setParam('TotalmRNAs',TotalmRNAs)
        self._setParam('mean_expression_threshold',mean_expression_threshold)
        self._setParam('num_PCA',num_PCA)
        self._setParam('mncluster_num',mncluster_num)
        #set sc3 parameters
        self._setParam('sc3matrix_file',sc3matrix_file)
        self._setParam('matrix_format',matrix_format)
        self._setParam('sc3cluster_num',sc3cluster_num)
        self._setParam('sc3ann',sc3ann)
        
    def _call(self,*args):
        # args[0]._getObj('FastqDump').getOutput('fastqOutput')
        pass
    
    def _build(self,):
        algorithm = self._getParam('algorithm')
        if algorithm=='Monocle':
            matrixdata = self._getParam('matrixdata')
            min_expression = self._getParam('min_expression')
            num_cells_expressed_threshold = self._getParam('num_cells_expressed_threshold')
            TotalmRNAs = self._getParam('TotalmRNAs')
            mean_expression_threshold = self._getParam('mean_expression_threshold')
            num_PCA = self._getParam('num_PCA')
            mncluster_num = self._getParam('mncluster_num')
            MonocleQC_result = MonocleQC(matrixdata = matrixdata, 
                                         outputpath = None,
                                         min_expression = min_expression,
                                         num_cells_expressed_threshold = num_cells_expressed_threshold,
                                         TotalmRNAs = TotalmRNAs, 
                                         mean_expression_threshold = mean_expression_threshold)
            MonocleDC_result = MonocleDC(outputpath = None,
                                         num_PCA = num_PCA,
                                         cluster_num = mncluster_num)(MonocleQC_result)
            rp = Report()
            rp.add('Monocle Quality Control Results',[MonocleQC_result])
            rp.add('Monocle Dimension Reduction And Density Peak Clustering Results',[MonocleDC_result])
            self._setObj('MonocleQC_result',MonocleQC_result)
            self._setObj('MonocleDC_result',MonocleDC_result)
        elif algorithm=='sc3':
            sc3matrix_file = self._getParam('sc3matrix_file')
            matrix_format = self._getParam('matrix_format')
            sc3cluster_num = self._getParam('sc3cluster_num')
            sc3ann = self._getParam('sc3ann')
            sce = SingleCellExperiment(matrix_file=sc3matrix_file,
                                       ann_file=sc3ann, 
                                       matrix_format='ORIGIN')
            sc3_cluster = SC3_Cluster(cluster_num = sc3cluster_num)(sce)
            rp = Report()
            rp.add('SingleCellExperiment Results',[sce])
            rp.add('SC3 Clustering Results',[sc3_cluster])
            self._setObj('sce',sce)
            self._setObj('sc3_cluster',sc3_cluster)
        else:
            print('This algorithm is unavailable!') 

    def _copy(self,):
        algorithm = self._getParam('algorithm')
        if algorithm=='Monocle':
            pass
        elif algorithm=='sc3':
            self._linkRecursive(self._getObj('sc3_cluster').getOutput('sc3Output_Result.xls'),
                                self.getFinalRsPath('sc3Output_Result'))
        else:
            print('This algorithm is unavailable!') 