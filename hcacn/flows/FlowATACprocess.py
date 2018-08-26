# -*- coding: utf-8 -*-
"""
@Time    : 2018/4/2 10:42
@Author  : Weizhang
@FileName: FlowATACprocess.py
"""

from ..core import Flow, Report, Configure, Schedule
from ..steps import SRAToFastq, AdapterRemoval, Bowtie2, SamToBam, BamSort, RmDuplicates, BamToBed
from ..steps import RmChrOrMergeAllSample, MergeToFrag, BedSort, PeakCalling, GenPeakWithFilter
from ..steps import LibComplexity, FragInPeak, CellFilter, CellExtracterBam, FragLenDistri


class FlowATACprocess(Flow):
    def __init__(self,
                 sraInput,  # for SRAToFatsq
                 refdir,  # for bowtie2
                 genome,  # for bowtie2
                 resultDir='./result',
                 blacklist=None,  # for GenPeakWithFilter
                 threads=1,
                 adapter1=None,  # for adapaterremoval
                 adapter2=None,
                 isdovetail=True,  # for bowtie2
                 isNoDiscordant=True,
                 isNoUnal=True,
                 isNoMixed=True,
                 X=2000,
                 memory='-Xmx5g',  # for RmDuplicates and LibComplexity
                 FSShift=4,  # for BamToBed
                 RSShift=-5,
                 savedchr=None,  # for RmChrOrMergeAllSample
                 formatPK='BED',  # for peakcalling
                 genomePK='hs',
                 nomodel=True,
                 nolambda=True,
                 keepdup='--keep-dup all',
                 callsummits=True,
                 overlapRate=0.2, # for GenPeakWithFilter
                 extendRange=250,
                 topPeak=50000,
                 libSizeCutOff=10000,  # CellFilter
                 fragInPeakCutOff=0.15,
                 ):
        super(FlowATACprocess, self).__init__(resultDir=resultDir,
                                              refdir=refdir,
                                              genome=genome,
                                              threads=threads)
        self._setParam('sraInput', sraInput)
        self._setParam('blacklist', blacklist)
        self._setParam('adapter1', adapter1)
        self._setParam('adapter2', adapter2)
        self._setParam('isdovetail', isdovetail)
        self._setParam('isNoDiscordant', isNoDiscordant)
        self._setParam('isNoUnal', isNoUnal)
        self._setParam('isNoMixed', isNoMixed)
        self._setParam('X', X)
        self._setParam('memory', memory)
        self._setParam('FSShift', FSShift)
        self._setParam('RSShift', RSShift)
        self._setParam('savedchr', savedchr)
        self._setParam('formatPK', formatPK)
        self._setParam('genomePK', genomePK)
        self._setParam('nomodel', nomodel)
        self._setParam('nolambda', nolambda)
        self._setParam('keepdup', keepdup)
        self._setParam('callsummits', callsummits)
        self._setParam('overlapRate', overlapRate)
        self._setParam('extendRange', extendRange)
        self._setParam('topPeak', topPeak)
        self._setParam('libSizeCutOff', libSizeCutOff)
        self._setParam('fragInPeakCutOff', fragInPeakCutOff)


    def _call(self,*args):
        pass


    def _build(self,):
        sraInput = self._getParam('sraInput')
        blacklist = self._getParam('blacklist')
        adapter1 = self._getParam('adapter1')
        adapter2 = self._getParam('adapter2')
        isdovetail = self._getParam('isdovetail')
        isNoDiscordant = self._getParam('isNoDiscordant')
        isNoUnal = self._getParam('isNoUnal')
        isNoMixed = self._getParam('isNoMixed')
        X = self._getParam('X')
        memory = self._getParam('memory')
        FSShift = self._getParam('FSShift')
        RSShift = self._getParam('RSShift')
        savedchr = self._getParam('savedchr')
        formatPK = self._getParam('formatPK')
        genomePK = self._getParam('genomePK')
        nomodel = self._getParam('nomodel')
        nolambda = self._getParam('nolambda')
        keepdup = self._getParam('keepdup')
        callsummits = self._getParam('callsummits')
        overlapRate = self._getParam('overlapRate')
        extendRange = self._getParam('extendRange')
        topPeak = self._getParam('topPeak')
        libSizeCutOff = self._getParam('libSizeCutOff')
        fragInPeakCutOff = self._getParam('fragInPeakCutOff')

        stf = SRAToFastq(sraInput=sraInput)
        adrm = AdapterRemoval(adapter1=adapter1, adapter2=adapter2)(stf)
        rs = Bowtie2(threads=self.threads, isdovetail=isdovetail, isNoDiscordant=isNoDiscordant,
                     isNoUnal=isNoUnal, isNoMixed=isNoMixed, X=X)(adrm)
        sb = SamToBam(threads=self.threads)(rs)
        bs = BamSort(threads=self.threads)(sb)
        rd = RmDuplicates(memory=memory)(bs)
        bb = BamToBed(FSShift=FSShift, RSShift=RSShift)(rd)
        rcmas = RmChrOrMergeAllSample(savedchr=savedchr)(bb)
        mtf2 = MergeToFrag()(rcmas)
        fld = FragLenDistri()(mtf2)
        bedsort = BedSort()(rcmas)
        peakc = PeakCalling(format=formatPK, genome=genomePK, nomodel=nomodel, nolambda=nolambda,
                            keepdup=keepdup, callsummits=callsummits)(bedsort)
        toppeak = GenPeakWithFilter(topPeak=topPeak, overlapRate=overlapRate,
                                    extendRange=extendRange, blacklist=blacklist)(peakc)
        lc = LibComplexity(memory=memory)(rd)
        fip = FragInPeak()(mtf2, toppeak)
        cf = CellFilter(libSizeCutOff=libSizeCutOff, fragInPeakCutOff=fragInPeakCutOff)(lc, fip)
        ceb = CellExtracterBam()(cf, rd)

        rp = Report()
        rp.add('Section for FASTQ Quality Control', [stf])
        rp.add('Section for Adapter Remove', [adrm])
        rp.add('Section for Mapping Result', [rs])
        rp.add('Section for Remove Chromatin', [rcmas])
        rp.add('Section for Fragment Length Distribution', [fld])
        rp.add('Section for Peak Calling', [peakc])
        rp.add('Section for Peak Filter', [toppeak])
        rp.add('Section for Cell Filter', [cf])

        self._setObj('SRAToFastq', stf)
        self._setObj('AdapterRemoval', adrm)
        self._setObj('Bowtie2', rs)
        self._setObj('SamToBam', sb)
        self._setObj('BamSort', bs)
        self._setObj('RmDuplicates', rd)
        self._setObj('BamToBed', bb)
        self._setObj('RmChrOrMergeAllSample', rcmas)
        self._setObj('MergeToFrag', mtf2)
        self._setObj('FragLenDistri', fld)
        self._setObj('BedSort', bedsort)
        self._setObj('PeakCalling', peakc)
        self._setObj('GenPeakWithFilter', toppeak)
        self._setObj('LibComplexity', lc)
        self._setObj('FragInPeak', fip)
        self._setObj('CellFilter', cf)
        self._setObj('CellExtracterBam', ceb)

    def _copy(self, ):
        self._linkRecursive(self._getObj('GenPeakWithFilter').getOutput('bedOutput'),
                            self.getFinalRsPath('TopPeak'))
        self._linkRecursive(self._getObj('CellExtracterBam').getOutput('bamOutput'),
                            self.getFinalRsPath('SavedBAM'))

