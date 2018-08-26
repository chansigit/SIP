# -*- coding: utf-8 -*-
"""
@Time    : 2018/4/2 14:48
@Author  : Weizhang
@FileName: test_FlowATAC.py
"""


from hcacn.core import Configure, Schedule
from hcacn.flows import FlowATACprocess, FlowATACAnalysis

chr_info = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

Configure.setIdentity('ATAC')

rs=FlowATACprocess(sraInput='/data8t_1/scATAC/GM12878',
                   refdir='/data8t_1/ref/atac/hg19_bowtie2',
                   genome='hg19',
                   threads=10,
                   resultDir='./result_ATAC',
                   blacklist='/data8t_1/ref/atac/others/hg19.blacklist.bed',
                   savedchr=chr_info)()

result=FlowATACAnalysis(refdir='/data8t_1/ref/atac/hg19_bowtie2',
                        genome='hg19',
                        threads=5,
                        resultDir='./result_ATAC')(rs)




