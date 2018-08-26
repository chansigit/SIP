# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
 """
from hcacn.core import Configure, Schedule
from hcacn.steps import Deseq2

Configure.setIdentity('zywang')
Deseq2(matrixdata = "./minidata/Deseq2/Deseq2testdata.txt", annotation = "./minidata/Deseq2/condition.csv", outputpath = None)
Schedule.run()