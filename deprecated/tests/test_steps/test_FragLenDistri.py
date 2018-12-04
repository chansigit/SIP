# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/30 19:39
@Author  : Weizhang
@FileName: test_FragLenDistri.py
"""

from hcacn.core import Configure,Schedule
from hcacn.steps import FragLenDistri

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

fld = FragLenDistri(bedInput='./minidata/atac/allbed')


Schedule.run()

