from hcacn.steps import Cellranger
from hcacn.core import Configure, Schedule


test = Cellranger(fastqInput = '/home/cfeng/data/test/', outputdir='test_cellranger', refile = '/home/cfeng/data/refdata-cellranger-hg19_and_mm10-1.2.0', expectcells=100)
Schedule.run()

print('')