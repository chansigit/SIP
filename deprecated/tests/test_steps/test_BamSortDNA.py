# -*- coding: utf-8 -*-

from hcacn.steps import BamSort
from hcacn.core import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=BamSort(bamInput='./minidata/atac/BamForTest')

Schedule.run()

