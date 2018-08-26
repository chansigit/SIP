# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 17:41
@Author  : Weizhang
@FileName: test_VarAndClustering.py
"""

from hcacn.steps import VarAndClustering
from hcacn.core import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

a=VarAndClustering(bamInput='./minidata/atac/bam_sorted_rmdup',
                   peakInput='./minidata/atac/others/top_peaks.bed',
                   threads=4,
                   genome='hg19')

Schedule.run()

