from hcacn.steps import Bamsort
from hcacn.core import Configure,Schedule

import os

Configure.setRefDir('./ref')
Configure.setGenome('hg19')

Bmst = Bamsort(bamInput='./hisat_fin.bam',
			   bamOutputDir='./')

Schedule.run()
