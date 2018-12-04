from hcacn.core import Configure,Schedule
from hcacn.steps import FastqDump
from hcacn.steps import Cuffmerge
# Configure.setRefDir('/home/zwei/ref')
# Configure.setGenome('hg19')


Cuffmerge(faInput1 = '/data/sqchen/hg19.fa',  
         gtfInput1 = '/data/sqchen/genome.gtf',  
         assembliesInput1 = '/data/sqchen/assemblies.txt',
         threads = 16,
         gtfOutputDir = '/data/sqchen')
Schedule.run()
 # docker run --rm -v /home/hca/Docker/Common_data:/data hca:py2 cuffmerge -g /data/sqchen/genome.gtf -s /data/sqchen/hg19.fa -o /data/sqchen -p 8 /data/sqchen/assemblies.txt
print('[done]')
