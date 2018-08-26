# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
 """
from hcacn.core import Configure, Schedule
from hcacn.steps import matrixRw
Configure.setIdentity('zywang')
matrixRw(matrixdata="./minidata/matrixRw/matrix.txt", outputpath=None)
Schedule.run()