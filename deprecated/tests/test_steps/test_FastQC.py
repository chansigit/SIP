
# coding: utf-8
from hcacn.core import Step,Configure,Schedule
from hcacn.steps import FastQC
Configure.setIdentity("yinqijin")

Configure.enableDocker(True)
# Folder Test
# fastqc = FastQC('./minidata/test_fastqc/','fastq',)
# FileTest
#fastqc = FastQC('./minidata/smartseq/fastq/','fastq',)
fastqc = FastQC('./minidata/smartseq/fastq/',)
Schedule.run()
