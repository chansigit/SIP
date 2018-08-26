from hcacn.steps import TrimPolyA
from hcacn.core import Configure, Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

ta = TrimPolyA(bamInput = './step_00_TrimAdapter/unaligned_tagged_trimmed_smart.bam', misMatches = 0, numBases = 6)

Schedule.run()
