from hcacn.steps import DetectError
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

de = DetectError(bamInput='./step_00_TagGene/star_gene_exon_tagged.bam',
                 numCells=100, primerSeqence='GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT')
Schedule.run()
