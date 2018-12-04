# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 11:28
@Author  : Weizhang
@FileName: test_GenPeakWithFilter.py
"""

from hcacn.core import Configure,Schedule
from hcacn.steps import GenPeakWithFilter

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

a=GenPeakWithFilter(summitInput='./minidata/atac/others/output_summits.bed',
                    blacklist='./minidata/atac/others/consensusBlacklist.bed',
                    topPeak=500)

Schedule.run()

