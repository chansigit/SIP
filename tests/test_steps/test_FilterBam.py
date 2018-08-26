from hcacn.steps import FilterBam
from hcacn.core import Configure, Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

fb = FilterBam(bamInput = './step_00_TagBarcode/unalign_tagged_XM.bam', tagReject = 'XQ')
Schedule.run()
