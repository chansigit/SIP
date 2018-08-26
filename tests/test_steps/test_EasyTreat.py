from hcacn.steps import EasyTreat
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

pt = EasyTreat(dgeInput='./step_13_DigitalExpression/out_gene_exon_tagged.dge.txt.gz')
Schedule.run()
