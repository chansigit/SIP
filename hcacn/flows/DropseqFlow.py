# -*- coding: utf-8 -*-

from ..core import Flow, Report, Configure, Schedule
from ..steps import FastqToBam, BamMerge, TagBarcode, FilterBam, TrimAdapter, TrimPolyA, BamToFastq, StarAlign, SortBam, MergeBamAlign, TagGene, DetectError, DigitalExpression

import os

class DropseqFlow(Flow):
    def __init__(self,
                 fastqInput1,
                 fastqInput2,
                 refInput = None,
                 refdir= None,
                 genome = None,
                 resultDir = './result',
                 cellBarcodeStart = 1,
                 cellBarcodeEnd = 16,
                 molecularBarcodeStart = 17,
                 molecularBarcodeEnd = 26,
                 cellBarcodeRead = 1,
                 cellNumBaseBelowQuality = 1,
                 molecularBarcodeRead = 1,
                 molecularNumBaseBelowQuality = 1,
                 cellBaseQuality = 10,
                 molecularBaseQuality = 10,
                 cellDiscardRead = False,
                 molecularDiscardRead = True,
                 adaperMisMatches = 0,
                 adaperNumBases = 5,
                 polyAMisMatches = 0,
                 polyANumBases = 6,
                 adapterSeq = 'AAGCAGTGGTATCAACGCAGAGTACATGGG',
                 starThreads = 16,
                 starOutFileNamePrefix = 'star',
                 secondAlign = False,
                 pairedRun = False,
                 expectedCellNum = 10000):
        super(DropseqFlow, self).__init__(resultDir = resultDir,
                                          refdir = refdir,
                                          genome = genome,
                                          threads = starThreads)
        self._setParam('fastqInput1', fastqInput1)
        self._setParam('fastqInput2', fastqInput2)
        self._setParam('adapterSeq', adapterSeq)
        self._setParam('cellBarcodeStart', cellBarcodeStart)
        self._setParam('cellBarcodeEnd', cellBarcodeEnd)
        self._setParam('molecularBarcodeStart', molecularBarcodeStart)
        self._setParam('molecularBarcodeEnd', molecularBarcodeEnd)
        self._setParam('starThreads', starThreads)
        self._setParam('starOutFileNamePrefix', starOutFileNamePrefix)
        self._setParam('expectedCellNum', expectedCellNum)
        self._setParam('refInput', refInput)
        self._setParam('cellBarcodeRead', cellBarcodeRead)
        self._setParam('cellNumBaseBelowQuality', cellNumBaseBelowQuality)
        self._setParam('molecularBarcodeRead', molecularBarcodeRead)
        self._setParam('molecularNumBaseBelowQuality', molecularNumBaseBelowQuality)
        self._setParam('cellDiscardRead', cellDiscardRead)
        self._setParam('molecularDiscardRead', molecularDiscardRead)
        self._setParam('cellBaseQuality', cellBaseQuality)
        self._setParam('molecularBaseQuality', molecularBaseQuality)
        self._setParam('adaperMisMatches', adaperMisMatches)
        self._setParam('adaperNumBases', adaperNumBases)
        self._setParam('polyAMisMatches', polyAMisMatches)
        self._setParam('polyANumBases', polyANumBases)
        self._setParam('secondAlign', secondAlign)
        self._setParam('pairedRun', pairedRun)
        print('DropseqFlow init')

    def _call(self, *args):
        pass

    def _build(self, ):
        fastqInput1 = self._getParam('fastqInput1')
        fastqInput2 = self._getParam('fastqInput2')
        adapterSeq = self._getParam('adapterSeq')
        cellBarcodeStart = self._getParam('cellBarcodeStart')
        cellBarcodeEnd = self._getParam('cellBarcodeEnd')
        molecularBarcodeStart = self._getParam('molecularBarcodeStart')
        molecularBarcodeEnd = self._getParam('molecularBarcodeEnd')
        starThreads = self._getParam('starThreads')
        starOutFileNamePrefix = self._getParam('starOutFileNamePrefix')
        expectedCellNum = self._getParam('expectedCellNum')
        refInput = self._getParam('refInput')
        cellBarcodeRead = self._getParam('cellBarcodeRead')
        cellNumBaseBelowQuality = self._getParam('cellNumBaseBelowQuality')
        molecularBarcodeRead = self._getParam('cellBarcodeRead')
        molecularNumBaseBelowQuality = self._getParam('molecularNumBaseBelowQuality')
        cellDiscardRead = self._getParam('cellDiscardRead')
        molecularDiscardRead = self._getParam('molecularDiscardRead')
        cellBaseQuality = self._getParam('cellBaseQuality')
        molecularBaseQuality = self._getParam('molecularBaseQuality')
        adaperMisMatches = self._getParam('adaperMisMatches')
        adaperNumBases = self._getParam('adaperNumBases')
        polyAMisMatches = self._getParam('polyAMisMatches')
        polyANumBases = self._getParam('polyANumBases')
        secondAlign = self._getParam('secondAlign')
        pairedRun = self._getParam('pairedRun')
        
        if refInput is None:
            starRef = None
            mergeRef = None
            gtfRef = None
        else:
            starRef = refInput + '/star'
            mergeRef = refInput + '/fasta'
            gtfRef = refInput + '/genes/genes.gtf'

        f2b = FastqToBam(fastqInput1 = fastqInput1, fastqInput2 = fastqInput2)
        bm = BamMerge()(f2b)
        tbc = TagBarcode(baseStart = cellBarcodeStart, baseEnd = cellBarcodeEnd,
                        barcodeRead = cellBarcodeRead, discardRead = cellDiscardRead,
                        tagName = 'XC', numBaseBelowQuality = cellNumBaseBelowQuality,
                        baseQuality = cellBaseQuality)(bm)
        tbm = TagBarcode(baseStart = molecularBarcodeStart, baseEnd = molecularBarcodeEnd,
                        barcodeRead = molecularBarcodeRead, discardRead = molecularDiscardRead,
                        tagName = 'XM', numBaseBelowQuality = molecularNumBaseBelowQuality,
                        baseQuality = cellBaseQuality)(tbc)
        fb = FilterBam(tagReject = 'XQ')(tbm)
        ta = TrimAdapter(adapterSeq = adapterSeq, misMatches = adaperMisMatches, numBases = adaperNumBases)(fb)
        tp = TrimPolyA(misMatches = polyAMisMatches, numBases = polyANumBases)(ta)
        b2f = BamToFastq()(tp)
        sa = StarAlign(outFileNamePrefix = starOutFileNamePrefix, genomeDir = starRef, threads = starThreads)(b2f)
        sb = SortBam(sortOrder = 'queryname')(sa) 
        mba = MergeBamAlign(refInputDir = mergeRef, secondAlign = secondAlign, pairedRun = pairedRun)(tp, sb)
        tg = TagGene(gtfInput = gtfRef, tag='GE')(mba)
        de = DetectError(numCells = expectedCellNum, primerSeqence = adapterSeq)(tg)
        dge = DigitalExpression(numCells = expectedCellNum)(de)
        
        rp = Report()
        rp.add('Section for FastqInput', [f2b])
        rp.add('Section for Tag Barcode', [tbc, tbm])
        rp.add('Section for Trim Reads', [ta, tp])
        rp.add('Section for STAR Alignment', [sa])
        rp.add('Section for DetectError', [de])
        rp.add('Section for DigitalExpression', [dge])

        self._setObj('FastqToBam', f2b)
        self._setObj('BamMerge', bm)
        self._setObj('TagCellBarcode', tbc)
        self._setObj('TagMolecularBarcode', tbm)
        self._setObj('FilterBam', fb)
        self._setObj('TrimAdapter', ta)
        self._setObj('TrimPolyA', tp)
        self._setObj('BamToFastq', b2f)
        self._setObj('StarAlign', sa)
        self._setObj('SortBam', sb)
        self._setObj('MergeBamAlign', mba)
        self._setObj('TagGene', tg)
        self._setObj('DetectError', de)
        self._setObj('DigitalExpression', dge)
        self._setObj('Report', rp)

        print('DropseqFlow builted')

    def _copy(self, ):
        self._linkRecursive(self._getObj('StarAlign').getOutput('logFinalOutput'),
                            self.getFinalRsDir())
        self._linkRecursive(self._getObj('StarAlign').getOutput('bamOutput'),
                            self.getFinalRsDir())
        self._linkRecursive(self._getObj('DigitalExpression').getOutput('dgeOutput'),
                            self.getFinalRsDir())
        self._linkRecursive(self._getObj('DigitalExpression').getOutput('sumOutput'),
                            self.getFinalRsDir())

