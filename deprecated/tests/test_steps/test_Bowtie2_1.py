from hcacn.steps import Bowtie2
from hcacn.steps import AdapterRemoval
from hcacn.core import Configure,Schedule

import os

Configure.setRefDir('/home/wzhang/genome/hg19_bowtie2/')
Configure.setGenome('hg19')
#Configure.setIdentity('zwei')
Configure.enableDocker(False)


rs=Bowtie2(fastqInput1='./minidata/atac/end1',
           fastqInput2='./minidata/atac/end2')

Schedule.run()


