

from ..core import Step,Configure
import os

class Cuffnorm(Step):
	def __init__(self,
				 gtfInput = None,
				 cxbInput = None,
				 outputDir = None,
				 markerInput = None,

				 threads = None,
				 cmdParam = None,
				 **kwargs
				 ):
		super(Step, self).__init__(cmdParam,**kwargs)

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
		gtfInput = self.getParamIO('gtfInput')
		cxbInput = self.getParamIO('cxbInput')
		markerInput = self.getParamIO('markerInput')
		outputDir = self.getParamIO('outputDir')
		if cxbInput is not None:
			for i,item in enumerate(cxbInput.split(' ')):
				self.setInput('cxbInput_%d'%i,item)


		if gtfInput is not None:
			self.setInput('gtfInput',gtfInput)
			
		if outputDir is None:
			self.setParamIO('outputDir',Configure.getTmpDir())

		self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))
		outputDir = self.getParamIO('outputDir')
		# self.setParamIO('cxbInput',cxbInput)


		isoforms_fpkm_table=list()
		genes_fpkm_table=list()
		cds_fpkm_table=list()
		tss_groups_fpkm_table=list()	
		isoforms_count_table=list()
		genes_count_table=list()
		cds_count_table=list()
		tss_groups_count_table=list()

		cds_attr_table=list()
		genes_attr_table=list()
		isoforms_attr_table=list()
		tss_groups_attr_table=list()
		run_info=list()
		samples_table=list()

		isoforms_fpkm_table.append(os.path.join(outputDir,'cuffnorm_'+str(0),'isoforms.fpkm_table'))
		genes_fpkm_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'genes.fpkm_table'))
		cds_fpkm_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'cds.fpkm_table'))
		tss_groups_fpkm_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'tss_groups.fpkm_table'))
		isoforms_count_table.append(os.path.join(outputDir,'cuffnorm_'+str(0),'isoforms.count_table'))
		genes_count_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'genes.count_table'))
		cds_count_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'cds.count_table'))
		tss_groups_count_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'tss_groups.count_table'))
		
		isoforms_attr_table.append(os.path.join(outputDir,'cuffnorm_'+str(0),'isoforms.attr_table'))
		genes_attr_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'genes.attr_table'))
		cds_attr_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'cds.attr_table'))
		tss_groups_attr_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'tss_groups.attr_table'))
		run_info.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'run.info'))
		samples_table.append(os.path.join(outputDir, 'cuffnorm_'+str(0),'samples.table'))

		self.setOutput('isoforms_fpkm_table',isoforms_fpkm_table)
		self.setOutput('genes_fpkm_table',genes_fpkm_table)
		self.setOutput('cds_fpkm_table',cds_fpkm_table)
		self.setOutput('tss_groups_fpkm_table',tss_groups_fpkm_table)
		self.setOutput('isoforms_count_table',isoforms_count_table)
		self.setOutput('genes_count_table',genes_count_table)
		self.setOutput('cds_count_table',cds_count_table)
		self.setOutput('tss_groups_count_table',tss_groups_count_table)

		self.setOutput('isoforms_attr_table',isoforms_attr_table)
		self.setOutput('genes_attr_table',genes_attr_table)
		self.setOutput('cds_attr_table',cds_attr_table)
		self.setOutput('tss_groups_attr_table',tss_groups_attr_table)
		self.setOutput('run_info',run_info)
		self.setOutput('samples_table',samples_table)

		self._setInputSize(1)
	def call(self, *args):
		cxbUpstream = args[0]
		gtfUpstream = args[1]
		cxb = cxbUpstream.getOutput('abundances_cxb')

		#print('===================')
		#print(cxb)		
		marker = ''
		for i in range(len(cxb)):
			#cxb[i] = self.convertToRealPath(cxb[i])
			marker = marker + cxb[i].split('/')[-2] + ','
		cxb = ' '.join(cxb)
		marker = marker[:-1]
		self.setParamIO('cxbInput',cxb)
		self.setParamIO('markerInput',marker)
		self.setParamIO('gtfInput',gtfUpstream.getOutput('merged_gtf'))

	def _singleRun(self,i):
		gtfInput = self.getParamIO('gtfInput')
		cxbInput = self.getParamIO('cxbInput')
		markerInput = self.getParamIO('markerInput')
		outputDir = self.getParamIO('outputDir')
		cxbFinPath = []
		for i,item in enumerate(cxbInput.split(' ')):
			cxbFinPath.append(self.getInput('cxbInput_'+str(i)))
		cxbFin = ' '.join(cxbFinPath)
		cmdline = [
				'cuffnorm',
				'-o',os.path.join(outputDir,'cuffnorm_'+str(0)),
				'-p',str(self.getParam('threads')),
				'-L',markerInput,
				gtfInput[0],
				cxbFin
				]

		result = self.callCmdline('V1', cmdline)
		f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
		f.write(result.stdout)
		f.close()

	def getMarkdownEN(self,):
		mdtext = """
### cuffnorm Result
The cuffnorm result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
Total map reads means that total number of reads mapped to genome
""".format(mapRs = self.getOutput('stdOutput'))

		return mdtext
