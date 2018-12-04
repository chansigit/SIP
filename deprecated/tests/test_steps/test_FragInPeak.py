# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/27 16:01
@Author  : Weizhang
@FileName: test_FragInPeak.py
"""


from hcacn.steps import FragInPeak
from hcacn.core import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=FragInPeak(fragInput='./minidata/atac/allbed',
                peakInput='./minidata/atac/others/top_peaks.bed')

Schedule.run()


