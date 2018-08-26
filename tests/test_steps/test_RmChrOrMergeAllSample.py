# -*- coding: utf-8 -*-

from hcacn.steps import RmChrOrMergeAllSample
from hcacn.core import Configure,Schedule

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

chr_info = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']


test=RmChrOrMergeAllSample(bedInput='./minidata/atac/BedForTest',
                           savedchr=chr_info)

Schedule.run()