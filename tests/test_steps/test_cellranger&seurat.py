from hcacn.steps import Cellranger
from hcacn.steps import Seuratpreprocessing
from hcacn.steps import Seuratrun
from hcacn.core import Configure, Schedule

Configure.setIdentity('fengchen')

test = Cellranger(fastqInput = '/home/hca/fengchen/data/10Xdata/fastqs/',  refile = '/home/hca/fengchen/data/10Xdata/refdata-cellranger-hg19_and_mm10-1.2.0', expectcells=100)
test2 = Seuratpreprocessing(rscript = '/home/hca/fengchen/zuoye/Seuratpreprocessing.R')(test)
test3 = Seuratrun(rscript='/home/hca/fengchen/zuoye/Seuratrun.R')(test2)

Schedule.run()

print('')