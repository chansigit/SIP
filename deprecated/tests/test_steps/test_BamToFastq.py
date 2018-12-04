from hcacn.steps import BamToFastq
from hcacn.core import Configure, Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

b2f = BamToFastq(bamInput = './step_00_TrimPolyA/unaligned_mc_tagged_polyA_filtered.bam')

Schedule.run()
