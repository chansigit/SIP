# -*- coding: utf-8 -*-
from hcacn.steps import SRAToFastq
from hcacn.core import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=SRAToFastq(sraInput='./minidata/atac/SraForTest')

Schedule.run()



