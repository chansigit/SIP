import sys

sys.path.append("/data8t_1/chenshengquan/zuoye")

from hcacn.core import Configure,Schedule
from hcacn.flows import FlowSmartseq

Configure.setIdentity('sqchen_f6')
### cuffquant
# rs=FlowSmartseq(sraInput='/data8t_1/chenshengquan/minidata/test_sra',alignMethod='hisat2', refdir='/data8t_1/ref/smartseq',genome='hg19',threads=16,resultDir='./resultFlowSmartseq_1')()
# rs=FlowSmartseq(sraInput='/data8t_1/chenshengquan/minidata/test_sra', alignMethod='tophat2', refdir='/data8t_1/ref/smartseq',genome='hg19',threads=16,resultDir='./resultFlowSmartseq_tophat2')()
### fastq input
# rs=FlowSmartseq(fastqInput1='/data8t_1/chenshengquan/minidata/test_fastq1', fastqInput2='/data8t_1/chenshengquan/minidata/test_fastq2', alignMethod='hisat2', refdir='/data8t_1/ref/smartseq',genome='hg19',threads=16,resultDir='./resultFlowSmartseq_Hisat2_fastq')()

### HTSeq
# rs=FlowSmartseq(sraInput='/data8t_1/chenshengquan/minidata/test_sra',alignMethod='hisat2', refdir='/data8t_1/ref/smartseq',genome='hg19',threads=16,resultType='htseq', resultDir='./resultFlowSmartseq_Hisat2_htseq')()
### fastq input
rs=FlowSmartseq(fastqInput1='/data8t_1/chenshengquan/minidata/test_fastq1', fastqInput2='/data8t_1/chenshengquan/minidata/test_fastq2', alignMethod='hisat2', refdir='/data8t_1/ref/smartseq',genome='hg19',threads=16,resultType='htseq',resultDir='./resultFlowSmartseq_Hisat2_fastq_htseq')()