from hcacn.steps import BamToBed
from hcacn.core import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

test=BamToBed(bamInput='./minidata/atac/BamForTest')

Schedule.run()

