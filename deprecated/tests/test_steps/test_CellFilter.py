# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/29 11:08
@Author  : Weizhang
@FileName: test_CellFilter.py
"""

from hcacn.steps import CellFilter
from hcacn.core import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')


cf = CellFilter(libCpxInput='./minidata/atac/libcpx',
                fragInPeakInput='./minidata/atac/others/MergedAllFiles_summits_filterd_FragInPeak.txt')

Schedule.run()

