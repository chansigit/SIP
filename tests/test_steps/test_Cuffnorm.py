from hcacn.steps import Cuffnorm
from hcacn.core import Configure,Schedule

import os

Configure.setRefDir('ref')
Configure.setGenome('hg19')

cfnm = Cuffnorm(outputdir='./',
				markerInput='Cuffquant.0.suffix',
				gtfInput='./minidata/genome.gtf',
				cxbInput='/data/songshaoming/qat_results/Cuffquant.0.suffix/abundances.cxb'
				)

Schedule.run()