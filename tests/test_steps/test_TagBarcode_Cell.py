from hcacn.steps import TagBarcode
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

bm = TagBarcode(bamInput='./step_00_BamMerge/unaligned_data.bam', baseStart = 1, baseEnd = 16, baseQuality = 10,
                barcodeRead = 1, discardRead = False, tagName = 'XC', numBaseBelowQuality = 1)
Schedule.run()
