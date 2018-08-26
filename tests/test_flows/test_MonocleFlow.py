# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
2018/3/30
"""
from hcacn.core import Configure,Schedule
from hcacn.flows import MonocleDimentionReductionAndClusterFlow

Configure.setIdentity('Monocleflowtest')

MonocleObj= MonocleDimentionReductionAndClusterFlow(matrixdata="./minidata/Monocle/out_gene_exon_tagged.dge.txt",
	                                                resultDir='./result')()