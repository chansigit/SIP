# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
 """
from hcacn.core import Configure,Schedule
from hcacn.steps import Monocle2QC,Monocle2Pseudo


Configure.setIdentity('yinqijin')
#Configure.enableDocker(False)


#Configure.enableDocker(False)
# there must be a column called groups in phenodata, if phenodata exist
MonocleQC_result = Monocle2QC(
	matrixdata='./minidata/downstream/monocle/exprs', 
	featuredata = './minidata/downstream/monocle/feature',
    phenodata = './minidata/downstream/monocle/pheno',
	outputpath="./step_my") 

# 
pseudo = Monocle2Pseudo()(MonocleQC_result)


Schedule.run()


