# -*- coding: utf-8 -*-

from ..core import Step,Configure
import os
import subprocess

class Cuffdiff(Step):

	def __init__(self,
			faInput = None,
			gtfInput = None,
			cxbInput = None,
			markerInput = None,
			outputDir = None,
			threads = None,
			cmdParam = None,
			**kwargs
			):
		super(Step, self).__init__(cmdParam,**kwargs)
		Configure.setRefSuffix('faInput','.fa',check=False)
		self.setParamIO('faInput',faInput)
		self.setParamIO('gtfInput',gtfInput)
		self.setParamIO('cxbInput',cxbInput)
		self.setParamIO('markerInput',markerInput)
		self.setParamIO('outputDir',outputDir)
		self.initIO()

		if threads is None:
			threads = Configure.getThreads()
		self.setParam('threads',threads)
		self._setUpstreamSize(2)


	def impInitIO(self,):
		faInput = self.getParamIO('faInput')
		gtfInput = self.getParamIO('gtfInput')
		cxbInput = self.getParamIO('cxbInput')
		markerInput = self.getParamIO('markerInput')
		outputDir = self.getParamIO('outputDir')
		if cxbInput is not None:
			for i,item in enumerate(cxbInput.split(' ')):
				self.setInput('cxbInput_%d'%i,item)

		if gtfInput is not None:
			self.setInput('gtfInput',gtfInput)

		if faInput is None:
			faInput = Configure.getConfig('faInput')
			self.setParamIO('faInput',faInput)
		self.setInput('faInput',faInput)
			
		if outputDir is None:
			self.setParamIO('outputDir',Configure.getTmpDir())
		self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))

		outputDir = self.getParamIO('outputDir')

		bias_params_info = list()
		cds_count_tracking=list()
		cds_diff=list()
		cds_exp_diff=list()
		cds_fpkm_tracking=list()
		cds_read_group_tracking=list()
		gene_exp_diff=list()
		genes_count_tracking=list()
		genes_fpkm_tracking=list()
		genes_read_group_tracking=list()
		isoform_exp_diff=list()
		isoforms_count_tracking=list()
		isoforms_fpkm_tracking=list()
		isoforms_read_group_tracking=list()
		promoters_diff=list()
		read_groups_info=list()
		run_info=list()
		splicing_diff=list()
		tss_group_exp_diff=list()
		tss_groups_count_tracking=list()
		tss_groups_fpkm_tracking=list()
		tss_groups_read_group_tracking=list()
		var_model_info=list()

		bias_params_info.append(os.path.join(outputDir,'cuffdiff_'+str(0),'bias_params.info'))
		cds_count_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'cds.count_tracking'))
		cds_diff.append(os.path.join(outputDir,'cuffdiff_'+str(0),'cds.diff'))
		cds_exp_diff.append(os.path.join(outputDir,'cuffdiff_'+str(0),'cds_exp.diff'))
		cds_fpkm_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'cds.fpkm_tracking'))
		cds_read_group_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'cds.read_group_tracking'))
		gene_exp_diff.append(os.path.join(outputDir,'cuffdiff_'+str(0),'gene_exp.diff'))
		genes_count_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'genes.count_tracking'))
		genes_fpkm_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'genes.fpkm_tracking'))
		genes_read_group_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'genes.read_group_tracking'))
		isoform_exp_diff.append(os.path.join(outputDir,'cuffdiff_'+str(0),'isoform_exp.diff'))
		isoforms_count_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'isoforms.count_tracking'))
		isoforms_fpkm_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'isoforms.fpkm_tracking'))
		isoforms_read_group_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'isoforms.read_group_tracking'))
		promoters_diff.append(os.path.join(outputDir,'cuffdiff_'+str(0),'promoters.diff'))
		read_groups_info.append(os.path.join(outputDir,'cuffdiff_'+str(0),'read_groups.info'))
		run_info.append(os.path.join(outputDir,'cuffdiff_'+str(0),'run.info'))
		splicing_diff.append(os.path.join(outputDir,'cuffdiff_'+str(0),'splicing.diff'))
		tss_group_exp_diff.append(os.path.join(outputDir,'cuffdiff_'+str(0),'tss_group_exp.diff'))
		tss_groups_count_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'tss_groups.count_tracking'))
		tss_groups_fpkm_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'tss_groups.fpkm_tracking'))
		tss_groups_read_group_tracking.append(os.path.join(outputDir,'cuffdiff_'+str(0),'tss_groups.read_group_tracking'))
		var_model_info.append(os.path.join(outputDir,'cuffdiff_'+str(0),'var_model.info'))

		
		self.setOutput('bias_params_info',bias_params_info)
		self.setOutput('cds_count_tracking',cds_count_tracking)
		self.setOutput('cds_diff',cds_diff)
		self.setOutput('cds_exp_diff',cds_exp_diff)
		self.setOutput('cds_fpkm_tracking',cds_fpkm_tracking)
		self.setOutput('cds_read_group_tracking',cds_read_group_tracking)
		self.setOutput('gene_exp_diff',gene_exp_diff)
		self.setOutput('genes_count_tracking',genes_count_tracking)
		self.setOutput('genes_fpkm_tracking',genes_fpkm_tracking)
		self.setOutput('genes_read_group_tracking',genes_read_group_tracking)
		self.setOutput('isoform_exp_diff',isoform_exp_diff)
		self.setOutput('isoforms_count_tracking',isoforms_count_tracking)
		self.setOutput('isoforms_fpkm_tracking',isoforms_fpkm_tracking)
		self.setOutput('isoforms_read_group_tracking',isoforms_read_group_tracking)
		self.setOutput('promoters_diff',promoters_diff)
		self.setOutput('read_groups_info',read_groups_info)
		self.setOutput('run_info',run_info)
		self.setOutput('splicing_diff',splicing_diff)
		self.setOutput('tss_group_exp_diff',tss_group_exp_diff)
		self.setOutput('tss_groups_count_tracking',tss_groups_count_tracking)
		self.setOutput('tss_groups_fpkm_tracking',tss_groups_fpkm_tracking)
		self.setOutput('tss_groups_read_group_tracking',tss_groups_read_group_tracking)
		self.setOutput('var_model_info',var_model_info)

		self._setInputSize(1)

	def call(self, *args):
		cxbUpstream = args[1]
		gtfUpstream = args[0]

		cxb = cxbUpstream.getOutput('abundances_cxb')
		marker = ''
		for i in range(len(cxb)):
			#cxb[i] = self.convertToRealPath(cxb[i])
			marker = marker + cxb[i].strip().split('/')[-2] + ','
		cxb = ' '.join(cxb)
		marker = marker[:-1]
		self.setParamIO('cxbInput',cxb)
		self.setParamIO('markerInput',marker)

		self.setParamIO('gtfInput',gtfUpstream.getOutput('merged_gtf'))

	def _singleRun(self,i):
		gtfInput = self.getParamIO('gtfInput')
		faInput = self.getParamIO('faInput')
		cxbInput = self.getParamIO('cxbInput')
		markerInput = self.getParamIO('markerInput')
		outputDir = self.getParamIO('outputDir')
		cxbFinPath = []
		for i,item in enumerate(cxbInput.split(' ')):
			cxbFinPath.append(self.getInput('cxbInput_'+str(i)))
		cxbFin = ' '.join(cxbFinPath)

		cmdline = [
				'cuffdiff',
				'-o',os.path.join(outputDir,'cuffdiff_' + str(0)),
				'-p',str(self.getParam('threads')),
				'-L',markerInput,
				'-b',faInput,
				gtfInput[0],
				cxbFin
				]
		result = self.callCmdline('V1', cmdline)
		f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
		f.write(result.stdout)
		f.close()

	def getMarkdownEN(self,):
		mdtext = """
### cuffdiff Result
The cuffdiff result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
Total map reads means that total number of reads mapped to genome
""".format(mapRs = self.getOutput('stdOutput'))

		return mdtext