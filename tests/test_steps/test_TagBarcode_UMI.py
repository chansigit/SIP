from hcacn.steps import TagBarcode
from hcacn.core import Configure,Schedule

import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

bm = TagBarcode(bamInput='./step_00_TagBarcode/unalign_tagged_XC.bam', baseStart = 17, baseEnd = 26, baseQuality = 10,
                barcodeRead = 1, discardRead = True, tagName = 'XM', numBaseBelowQuality = 1)
Schedule.run()
