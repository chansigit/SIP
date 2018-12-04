# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/29 15:01
@Author  : Weizhang
@FileName: test_LibComplexity.py
"""
from hcacn.steps import LibComplexity
from hcacn.core import Configure, Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

lcpx = LibComplexity(bamInput='./minidata/atac/BamForTest')

Schedule.run()