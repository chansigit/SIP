from hcacn.steps import SortBam
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

sb = SortBam(bamInput='./step_00_StarAlign/starAligned.out.sam', sortOrder = 'queryname')
Schedule.run()
