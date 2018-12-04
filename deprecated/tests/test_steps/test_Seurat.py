from hcacn.core import Configure,Schedule
from hcacn.steps import Seurat
import os

test = Seurat(outputdir = 'test_seurat', rscript = '/data/test_celranger/Seurat.R')
Schedule.run()

print('')
