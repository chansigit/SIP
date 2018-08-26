from hcacn.steps import DigitalExpression
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

de = DigitalExpression(bamInput='./step_00_DetectError/out_gene_exon_tagged.bam', dgeOutputDir='./minidata/dropseq/tmp/',
                 sumOutputDir='./minidata/dropseq/tmp', numCells=100)
Schedule.run()
