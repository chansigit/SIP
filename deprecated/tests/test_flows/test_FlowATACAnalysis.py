# -*- coding: utf-8 -*-
"""
@Time    : 2018/4/3 15:12
@Author  : Weizhang
@FileName: test_FlowATACAnalysis.py
"""

from hcacn.core import Configure, Schedule
from hcacn.flows import FlowATACAnalysis

Configure.setIdentity('ATAC')

result=FlowATACAnalysis(bamInput='./minidata/atac/selectedBam',
                        peakInput='./minidata/atac/others/top_peaks.bed',
                        refdir='/data8t_1/ref/atac/hg19_bowtie2',
                        genome='hg19',
                        threads=10,
                        resultDir='./result_ATAC')()

