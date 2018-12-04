# -*- coding: utf-8 -*-
from hcacn.steps import Bowtie2, AdapterRemoval
from hcacn.core  import Configure, Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

adrm = AdapterRemoval(fastqInput1='./minidata/atac/end1', fastqInput2='./minidata/atac/end2')

Schedule.run()

