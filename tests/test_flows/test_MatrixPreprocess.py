# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
2018/3/29
"""
from hcacn.core import Configure,Schedule
from hcacn.flows import MatrixPreprocess

Configure.setIdentity('MatrixPreprocessflowtest')

matObj= MatrixPreprocess(matrixdata="./minidata/matrixRW/matrix.txt", outputpath=None)()


