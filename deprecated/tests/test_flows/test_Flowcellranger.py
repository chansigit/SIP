from hcacn.core import Configure, Schedule
from hcacn.flows import FlowCellranger
Configure.setIdentity('Fengchen')

cellranger=FlowCellranger(fastqInput='/data8t_1/fengchen/10Xdata/fastqs',refdir='/data8t_1/ref/dropseq', genome='hg19', resultDir='./result')()
