# -*- coding: utf-8 -*-
from hcacn.core import Configure,Schedule
from hcacn.flows import FlowExample
Configure.setIdentity('weizheng1')

rs=FlowExample(fastqInput1='./minidata/atac/end1', fastqInput2='./minidata/atac/end2',refdir='/data8t_1/hca/ref/hg19_bowtie2',genome='hg19',threads=4,resultDir='result')()


