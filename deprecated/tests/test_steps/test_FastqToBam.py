from hcacn.steps import FastqToBam
from hcacn.core import Configure, Schedule

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

b2s = FastqToBam(fastqInput1='./minidata/dropseq/read1', fastqInput2='./minidata/dropseq/read2')
Schedule.run()

print('')
