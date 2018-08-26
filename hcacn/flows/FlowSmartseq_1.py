# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 11:12:29 2018

@author: sqchen
"""
from ..core import Flow, Report, Configure, Schedule
from ..steps import FastqDump, Hisat2, Tophat2, SamToBam, BamSort, Cufflinks, Cuffmerge, Cuffquant, Cuffnorm, Cuffdiff, FastQC, Stringtie
class FlowSmartseq_1(Flow):
    def __init__(self,
                 sraInput,
                 refdir, 
                 genome, 
                 resultDir='./result',                  
                 threads=None,):
        super(FlowSmartseq_1,self).__init__(resultDir=resultDir, 
                                         refdir=refdir, 
                                         genome=genome, 
                                         threads=threads)
        self._setParam('sraInput',sraInput)
        
    def _call(self,*args):
        # args[0]._getObj('FastqDump').getOutput('fastqOutput')
        pass
    
    def _build(self,):
        sraInput = self._getParam('sraInput')
        fastqdump = FastqDump(sraInput1=sraInput)
        hisat = Hisat2()(fastqdump)
        sam2bam = SamToBam()(hisat)
        bamsort = BamSort()(sam2bam)
        stringtie = Stringtie()(bamsort)
        cuffmerge = Cuffmerge()(stringtie)
        cuffquant = Cuffquant()(bamsort,cuffmerge)
        cuffdiff = Cuffdiff()(cuffmerge,cuffquant)
                
        rp = Report()
        rp.add('Section for FastqDump',[fastqdump])
        rp.add('Section for Hisat2',[hisat])
        rp.add('Section for SamToBam',[sam2bam])
        rp.add('Section for BamSort',[bamsort])
        rp.add('Section for Stringtie',[stringtie])
        rp.add('Section for Cuffmerge',[cuffmerge])
        rp.add('Section for Cuffquant',[cuffquant])
        rp.add('Section for Cuffdiff',[cuffdiff])
        
        self._setObj('FastqDump',fastqdump)
        self._setObj('Hisat2',hisat)
        self._setObj('SamToBam',sam2bam)
        self._setObj('BamSort',bamsort)
        self._setObj('Stringtie',stringtie)
        self._setObj('Cuffmerge',cuffmerge)
        self._setObj('Cuffquant',cuffquant)
        self._setObj('Cuffdiff',cuffdiff)
        self._setObj('Report',rp)

        
    def _copy(self,):
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('bias_params_info'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_count_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_exp_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_fpkm_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('cds_read_group_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('gene_exp_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('genes_count_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('genes_fpkm_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('genes_read_group_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('isoform_exp_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('isoforms_count_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('isoforms_fpkm_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('isoforms_read_group_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('promoters_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('read_groups_info'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('run_info'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('splicing_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('tss_group_exp_diff'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('tss_groups_count_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('tss_groups_fpkm_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('tss_groups_read_group_tracking'), self.getFinalRsPath('Cuffdiff_Res'))
        self._linkRecursive(self._getObj('Cuffdiff').getOutput('var_model_info'), self.getFinalRsPath('Cuffdiff_Res'))
