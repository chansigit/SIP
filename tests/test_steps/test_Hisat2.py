from hcacn.core import Schedule
from hcacn.steps import Hisat2



# hisat = Hisat2("../smartseq2_hisat/SRR1294845_1.fastq","../smartseq2_hisat/SRR1294845_2.fastq","../hisat/genome_tran",
#                "../smartseq2_hisat_result")

hisat = Hisat2("../smartseq2_hisat/SRR1295221_1.fastq","../smartseq2_hisat/SRR1295221_2.fastq","../hg19_index/genome",
               "../smartseq2_hisat_result")

hisat.outputs


hisat.params


Schedule.run()

