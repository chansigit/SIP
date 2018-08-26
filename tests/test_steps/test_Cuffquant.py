from hcacn.steps import Cuffquant
from hcacn.core import Configure,Schedule

import os

Configure.setRefDir('ref')
Configure.setGenome('hg19')

cflk = Cuffquant(bamInput='./minidata/accepted_hits.bam',
				gtfInput='./minidata/genome.gtf',
				outputDir='./'
				)

Schedule.run()