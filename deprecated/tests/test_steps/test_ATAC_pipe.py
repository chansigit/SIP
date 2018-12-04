# -*- coding: utf-8 -*-

from hcacn.core import Configure, Schedule
from hcacn.steps import *
#from hcacn.steps import Bowtie2, AapterRemoval, SamToBam, BamSort, RmDuplicates,  BamToBed
#from hcacn.steps import RmChrOrMergeAllSample,MergeToFrag,BedSort,SRAToFastq, PeakCalling
#from hcacn.steps  import GenPeakWithFilter, VarAndClustering,LibComplexity,FragInPeak,CellFilter

Configure.setRefDir('/data8t_1/hca/ref/hg19_bowtie2')
Configure.setGenome('hg19')
Configure.setIdentity('ATAC')

stf = SRAToFastq(sraInput='/data8t_1/scATAC/GM12878')

adrm = AdapterRemoval()(stf)

rs=Bowtie2(threads=6)(adrm)

sb=SamToBam(threads=6)(rs)

bs=BamSort(threads=6)(sb)

rd=RmDuplicates(memory='-Xmx4g')(bs)

bb=BamToBed()(rd)

# # remain useless chromatin
# mtf1=MergeToFrag()(bb)

chr_info = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
            'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

rcmas=RmChrOrMergeAllSample(savedchr=chr_info)(bb)

mtf2=MergeToFrag()(rcmas)

fld=FragLenDistri()(mtf2)

bedsort=BedSort()(rcmas)

peakc=PeakCalling()(bedsort)

toppeak=GenPeakWithFilter(topPeak=50000,
                          blacklist='./minidata/atac/others/consensusBlacklist.bed')(peakc)

lc=LibComplexity(memory='-Xmx4g')(rd)

fip=FragInPeak()(mtf2, toppeak)

cf=CellFilter()(lc, fip)

ceb=CellExtracterBam()(cf, rd)

result=VarAndClustering(threads=3)(ceb, toppeak)

Schedule.run()
