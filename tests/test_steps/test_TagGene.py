from hcacn.steps import TagGene
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

tg = TagGene(bamInput='./step_00_MergeBamAlign/merged.bam',
                gtfInput = '../ref/refdata-cellranger-hg19_and_mm10-2.1.0/genes/genes.gtf',
                tag='GE')
Schedule.run()
