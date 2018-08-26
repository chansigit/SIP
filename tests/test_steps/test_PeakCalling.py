# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/12 21:52
@Author  : Weizhang
@FileName: test_PeakCalling.py
"""

from hcacn.steps import PeakCalling
from hcacn.core import Configure, Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=PeakCalling(bedInput='./minidata/atac/others/MergedAllFiles.bed')

Schedule.run()


