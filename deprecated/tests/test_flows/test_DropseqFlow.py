# -*- coding: utf-8 -*-
from hcacn.core import Configure,Schedule
from hcacn.flows import DropseqFlow
Configure.setIdentity('lcy')

df=DropseqFlow(fastqInput1='./minidata/dropseq/read1', fastqInput2='./minidata/dropseq/read2',refInput='/data8t_1/ref/dropseq/hg19-and-mm10_refdata-cellranger-2.1.0', refdir='/data8t_1/ref/dropseq/', genome='hg19-and-mm10', expectedCellNum=20000)()

