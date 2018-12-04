from hcacn.steps import MergeBamAlign
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

mba = MergeBamAlign(unmappedBamInput='./step_00_TrimPolyA/unaligned_mc_tagged_polyA_filtered.bam', alignedBamInput='./step_00_SortBam/aligned.sorted.bam',
                    refInputDir='../ref/refdata-cellranger-hg19_and_mm10-2.1.0/fasta/',
                    secondAlign=False, pairedRun=False)
Schedule.run()
