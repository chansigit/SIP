from hcacn.steps import BamMerge
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

bm = BamMerge(bamInput=['./step_00_FastqToBam/sample.0.bam', './step_00_FastqToBam/sample.1.bam'])
Schedule.run()
