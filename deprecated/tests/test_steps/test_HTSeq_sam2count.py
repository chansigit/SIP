from hcacn.core import Configure,Schedule
from hcacn.steps import FastqDump
from hcacn.steps import HTSeq_sam2count
# Configure.setRefDir('/home/zwei/ref')
# Configure.setGenome('hg19')


HTSeq_sam2count(samInput1='./accepted_hits.sam', gtfInput1='./genome.gtf')
Schedule.run()

print('[done]')
