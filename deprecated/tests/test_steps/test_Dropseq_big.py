from hcacn.steps import FastqToBam
from hcacn.steps import BamMerge
from hcacn.steps import TagBarcode
from hcacn.steps import FilterBam
from hcacn.steps import TrimAdapter
from hcacn.steps import TrimPolyA
from hcacn.steps import BamToFastq
from hcacn.steps import StarAlign
from hcacn.steps import SortBam
from hcacn.steps import MergeBamAlign
from hcacn.steps import TagGene
from hcacn.steps import DetectError
from hcacn.steps import DigitalExpression
from hcacn.steps import EasyTreat
from hcacn.steps import MonocleQC
from hcacn.steps import Monocle_dimreduce_cluster
from hcacn.core import Configure,Schedule
import os

#Configure.enableDocker(False)
Configure.setIdentity('cyliu')

f2b = FastqToBam(fastqInput1='../data/nuclei_900_fastqs/read1', fastqInput2='../data/nuclei_900_fastqs/read2')
                 #fastqInput1 = './minidata/dropseq/read1', fastqInput2 = './minidata/dropseq/read2')
                 #fastqInput1 = '../data/hgmm_100/read1', fastqInput2 = '../data/hgmm_100/read2')
bm = BamMerge()(f2b)
tbc = TagBarcode(baseStart = 1, baseEnd = 16, baseQuality = 10,
                barcodeRead = 1, discardRead = False, tagName = 'XC', numBaseBelowQuality = 1)(bm)
tbm = TagBarcode(baseStart = 17, baseEnd = 26, baseQuality = 10,
                barcodeRead = 1, discardRead = True, tagName = 'XM', numBaseBelowQuality = 1)(tbc)
fb = FilterBam(tagReject = 'XQ')(tbm)
ta = TrimAdapter(adapterSeq = 'AAGCAGTGGTATCAACGCAGAGTACATGGG',
                  misMatches = 0, numBases = 5)(fb)
tp = TrimPolyA(misMatches = 0, numBases = 6)(ta)
b2f = BamToFastq()(tp)
sa = StarAlign(outFileNamePrefix='star', genomeDir = '../ref/refdata-cellranger-mm10-2.1.0/star/', threads=16)(b2f)
sb = SortBam(sortOrder = 'queryname')(sa)
mba = MergeBamAlign(refInputDir='../ref/refdata-cellranger-mm10-2.1.0/fasta/',
                secondAlign=False, pairedRun=False)(tp, sb)
tg = TagGene(gtfInput = '../ref/refdata-cellranger-mm10-2.1.0/genes/genes.gtf', tag='GE')(mba)
de = DetectError(numCells=900, primerSeqence='AAGCAGTGGTATCAACGCAGAGTACATGGG')(tg)
dge = DigitalExpression(numCells=900)(de)
#et = EasyTreat()(dge)
mq = MonocleQC()(dge)
mdc = Monocle_dimreduce_cluster(cluster_num=5)(mq) 
Schedule.run()
