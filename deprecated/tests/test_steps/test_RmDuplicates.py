# -*- coding: utf-8 -*-
from hcacn.steps import RmDuplicates
from hcacn.core import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=RmDuplicates(bamInput='./minidata/atac/BamForTest',
                  memory='-Xmx4g')

Schedule.run()

