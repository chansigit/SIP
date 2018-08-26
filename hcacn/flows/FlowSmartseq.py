# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 11:12:29 2018

@author: sqchen
"""
from ..core import Flow, Report, Configure, Schedule
from ..steps import FastqDump, Hisat2, Tophat2, SamToBam, BamSort, Cufflinks, Cuffmerge, Cuffquant, Cuffnorm, Cuffdiff, FastQC, Stringtie, HTSeq_sam2count
class FlowSmartseq(Flow):
    def __init__(self,
                 sraInput=None,
                 fastqInput1=None,
                 fastqInput2=None,
                 refdir=None, 
                 genome=None,
                 resultType='cuffquant',
                 resultDir='./result',
                 alignMethod='hisat2',
                 threads=None,):
        super(FlowSmartseq,self).__init__(resultDir=resultDir, 
                                         refdir=refdir, 
                                         genome=genome, 
                                         threads=threads)
        self._setParam('sraInput',sraInput)
        self._setParam('fastqInput1',fastqInput1)
        self._setParam('fastqInput2',fastqInput2)
        self._setParam('alignMethod',alignMethod)
        self._setParam('resultType',resultType)
        
    def _call(self,*args):
        # args[0]._getObj('FastqDump').getOutput('fastqOutput')
        pass
    
    def _build(self,):
        sraInput = self._getParam('sraInput')
        fastqInput1 = self._getParam('fastqInput1')
        fastqInput2 = self._getParam('fastqInput2')
        alignMethod = self._getParam('alignMethod')
        resultType = self._getParam('resultType')
        if resultType == 'cuffquant':
            if sraInput is not None:
                if alignMethod == 'tophat2' :
                    fastqdump = FastqDump(sraInput1=sraInput)
                    tophat = Tophat2()(fastqdump)
                    sam2bam = SamToBam()(tophat)
                    bamsort = BamSort()(sam2bam)
                    cufflinks = Cufflinks()(bamsort)
                    cuffmerge = Cuffmerge()(cufflinks)
                    cuffquant = Cuffquant()(bamsort,cuffmerge)
                            
                    rp = Report()
                    rp.add('Section for FastqDump',[fastqdump])
                    rp.add('Section for Tophat2',[tophat])
                    rp.add('Section for SamToBam',[sam2bam])
                    rp.add('Section for BamSort',[bamsort])
                    rp.add('Section for Cufflinks',[cufflinks])
                    rp.add('Section for Cuffmerge',[cuffmerge])
                    rp.add('Section for Cuffquant',[cuffquant])
                    
                    self._setObj('FastqDump',fastqdump)
                    self._setObj('Tophat2',tophat)
                    self._setObj('SamToBam',sam2bam)
                    self._setObj('BamSort',bamsort)
                    self._setObj('Cufflinks',cufflinks)
                    self._setObj('Cuffmerge',cuffmerge)
                    self._setObj('Cuffquant',cuffquant)
                    self._setObj('Report',rp)
                else:
                    fastqdump = FastqDump(sraInput1=sraInput)
                    hisat = Hisat2()(fastqdump)
                    sam2bam = SamToBam()(hisat)
                    bamsort = BamSort()(sam2bam)
                    stringtie = Stringtie()(bamsort)
                    cuffmerge = Cuffmerge()(stringtie)
                    cuffquant = Cuffquant()(bamsort,cuffmerge)
                    # cuffdiff = Cuffdiff()(cuffmerge,cuffquant)
                            
                    rp = Report()
                    rp.add('Section for FastqDump',[fastqdump])
                    rp.add('Section for Hisat2',[hisat])
                    rp.add('Section for SamToBam',[sam2bam])
                    rp.add('Section for BamSort',[bamsort])
                    rp.add('Section for Stringtie',[stringtie])
                    rp.add('Section for Cuffmerge',[cuffmerge])
                    rp.add('Section for Cuffquant',[cuffquant])
                    # rp.add('Section for Cuffdiff',[cuffdiff])
                    
                    self._setObj('FastqDump',fastqdump)
                    self._setObj('Hisat2',hisat)
                    self._setObj('SamToBam',sam2bam)
                    self._setObj('BamSort',bamsort)
                    self._setObj('Stringtie',stringtie)
                    self._setObj('Cuffmerge',cuffmerge)
                    self._setObj('Cuffquant',cuffquant)
                    # self._setObj('Cuffdiff',cuffdiff)
                    self._setObj('Report',rp)
            ### fastq input
            else:
                if alignMethod == 'tophat2' :
                    tophat = Tophat2(fastqInput1=fastqInput1, fastqInput2=fastqInput2)
                    sam2bam = SamToBam()(tophat)
                    bamsort = BamSort()(sam2bam)
                    cufflinks = Cufflinks()(bamsort)
                    cuffmerge = Cuffmerge()(cufflinks)
                    cuffquant = Cuffquant()(bamsort,cuffmerge)
                            
                    rp = Report()
                    rp.add('Section for Tophat2',[tophat])
                    rp.add('Section for SamToBam',[sam2bam])
                    rp.add('Section for BamSort',[bamsort])
                    rp.add('Section for Cufflinks',[cufflinks])
                    rp.add('Section for Cuffmerge',[cuffmerge])
                    rp.add('Section for Cuffquant',[cuffquant])
                    
                    self._setObj('Tophat2',tophat)
                    self._setObj('SamToBam',sam2bam)
                    self._setObj('BamSort',bamsort)
                    self._setObj('Cufflinks',cufflinks)
                    self._setObj('Cuffmerge',cuffmerge)
                    self._setObj('Cuffquant',cuffquant)
                    self._setObj('Report',rp)
                else:
                    hisat = Hisat2(fastqInput1=fastqInput1, fastqInput2=fastqInput2)
                    sam2bam = SamToBam()(hisat)
                    bamsort = BamSort()(sam2bam)
                    stringtie = Stringtie()(bamsort)
                    cuffmerge = Cuffmerge()(stringtie)
                    cuffquant = Cuffquant()(bamsort,cuffmerge)
                            
                    rp = Report()
                    rp.add('Section for Hisat2',[hisat])
                    rp.add('Section for SamToBam',[sam2bam])
                    rp.add('Section for BamSort',[bamsort])
                    rp.add('Section for Stringtie',[stringtie])
                    rp.add('Section for Cuffmerge',[cuffmerge])
                    rp.add('Section for Cuffquant',[cuffquant])
                    
                    self._setObj('Hisat2',hisat)
                    self._setObj('SamToBam',sam2bam)
                    self._setObj('BamSort',bamsort)
                    self._setObj('Stringtie',stringtie)
                    self._setObj('Cuffmerge',cuffmerge)
                    self._setObj('Cuffquant',cuffquant)
                    self._setObj('Report',rp)

        ########## resultType == 'HTSeq'
        else:
            if sraInput is not None:
                if alignMethod == 'tophat2' :
                    fastqdump = FastqDump(sraInput1=sraInput)
                    tophat = Tophat2()(fastqdump)
                    htseq = HTSeq_sam2count()(tophat)
                            
                    rp = Report()
                    rp.add('Section for FastqDump',[fastqdump])
                    rp.add('Section for Tophat2',[tophat])
                    rp.add('Section for HTSeq_sam2count',[htseq])
                    
                    self._setObj('FastqDump',fastqdump)
                    self._setObj('Tophat2',tophat)
                    self._setObj('HTSeq_sam2count',htseq)
                    self._setObj('Report',rp)
                else:
                    fastqdump = FastqDump(sraInput1=sraInput)
                    hisat = Hisat2()(fastqdump)
                    htseq = HTSeq_sam2count()(hisat)
                            
                    rp = Report()
                    rp.add('Section for FastqDump',[fastqdump])
                    rp.add('Section for Hisat2',[hisat])
                    rp.add('Section for HTSeq_sam2count',[htseq])
                    
                    self._setObj('FastqDump',fastqdump)
                    self._setObj('Hisat2',hisat)
                    self._setObj('HTSeq_sam2count',htseq)
                    self._setObj('Report',rp)
            ### fastq input
            else:
                if alignMethod == 'tophat2' :
                    tophat = Tophat2(fastqInput1=fastqInput1, fastqInput2=fastqInput2)
                    htseq = HTSeq_sam2count()(tophat)
                            
                    rp = Report()
                    rp.add('Section for Tophat2',[tophat])
                    rp.add('Section for HTSeq_sam2count',[htseq])
                    
                    self._setObj('Tophat2',tophat)
                    self._setObj('HTSeq_sam2count',htseq)
                    self._setObj('Report',rp)
                else:
                    hisat = Hisat2(fastqInput1=fastqInput1, fastqInput2=fastqInput2)
                    htseq = HTSeq_sam2count()(hisat)
                            
                    rp = Report()
                    rp.add('Section for Hisat2',[hisat])
                    rp.add('Section for HTSeq_sam2count',[htseq])
                    
                    self._setObj('Hisat2',hisat)
                    self._setObj('HTSeq_sam2count',htseq)
                    self._setObj('Report',rp)

        
    def _copy(self,):
        resultType = self._getParam('resultType')
        if resultType == 'cuffquant':
            self._linkRecursive(self._getObj('Cuffquant').getOutput('abundances_cxb'), self.getFinalRsPath('Cuffquant_Res'))
        else:
            self._linkRecursive(self._getObj('HTSeq_sam2count').getOutput('countOutput'), self.getFinalRsPath('HTSeq_sam2count_Res'))
            
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('bias_params_info'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_count_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_exp_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_fpkm_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_read_group_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('gene_exp_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('genes_count_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('genes_fpkm_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('genes_read_group_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('isoform_exp_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('isoforms_count_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('isoforms_fpkm_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('isoforms_read_group_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('promoters_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('read_groups_info'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('run_info'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('splicing_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('tss_group_exp_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('tss_groups_count_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('tss_groups_fpkm_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('tss_groups_read_group_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        # self._linkRecursive(self._getObj('Cuffdiff').getOutput('var_model_info'), self.getFinalRsPath('Cuffdiff_Res'))
