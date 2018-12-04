from hcacn.steps import TrimAdapter
from hcacn.core import Configure, Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

#ta = TrimAdapter(bamInput = './step_00_FilterBam/unalign_tagged_filterd.bam', adapterSeq = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT', misMatches = 0, numBases = 5)
ta = TrimAdapter(bamInput = './step_00_FilterBam/unalign_tagged_filterd.bam', adapterSeq = 'AAGCAGTGGTATCAACGCAGAGTACATGGG', misMatches = 0, numBases = 5)
Schedule.run()
