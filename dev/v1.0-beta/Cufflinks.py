

from ..core import Step,Configure
import os

class Cufflinks(Step):
	Configure.setRefSuffix('gtfInput','.gtf',check=False)

	def __init__(self,
				 bamInput = None,
				 gtfInput = None,
				 outputDir = None,
				 threads = None,
				 ismultiReadCorrect = None,
				 isupperQuartileForm = None,
				 istotalHitsNorm = True,
				 fragLenMean = 200,
				 fragLenStdDev = 80,
				 cmdParam = None,
				 **kwargs
				):
		super(Step, self).__init__(cmdParam,**kwargs)
		self.setParamIO('bamInput',bamInput)
		self.setParamIO('gtfInput',gtfInput)
		self.setParamIO('outputDir',outputDir)
		#self.setParamIO('fragBiasCorrectInput',fragBiasCorrectInput)
		self.initIO()

		self.setParam('ismultiReadCorrect',ismultiReadCorrect)
		self.setParam('fragLenMean',fragLenMean)
		self.setParam('fragLenStdDev',fragLenStdDev)
		self.setParam('isupperQuartileForm',isupperQuartileForm)
		self.setParam('istotalHitsNorm',istotalHitsNorm)

		if threads is None:
			threads = Configure.getThreads()
		self.setParam('threads',threads)

	def impInitIO(self,):
		bamInput = self.getParamIO('bamInput')
		gtfInput = self.getParamIO('gtfInput')
		outputDir = self.getParamIO('outputDir')
		#fragBiasCorrectInput = self.getParamIO('fragBiasCorrectInput')
		if outputDir is None:
			self.setParamIO('outputDir',Configure.getTmpDir())

		self.setInputDirOrFile('bamInput',bamInput)

		if gtfInput is None:
			gtfInput = Configure.getConfig('gtfInput')
			self.setParamIO('gtfInput',gtfInput)
		self.setInput('gtfInput',gtfInput)

		self.setOutput('assembliesOutput',os.path.join(Configure.getTmpDir(), 'assemblies.txt'))
		self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))
		
		if bamInput is not None:
			self._setInputSize(len(self.getInputList('bamInput')))
			genes_fpkm_tracking=list()
			isoforms_fpkm_tracking=list()
			skipped_gtf=list()
			transcripts_gtf=list()
			for i in range(len(self.getInputList('bamInput'))):
				genes_fpkm_tracking.append(os.path.join(outputDir, 'cufflinks_'+str(i),'genes.fpkm_tracking'))
				isoforms_fpkm_tracking.append(os.path.join(outputDir, 'cufflinks_'+str(i),'isoforms.fpkm_tracking'))
				skipped_gtf.append(os.path.join(outputDir, 'cufflinks_'+str(i),'skipped.gtf'))
				transcripts_gtf.append(os.path.join(outputDir, 'cufflinks_'+str(i),'transcripts.gtf'))
			self.setOutput('genes_fpkm_tracking',genes_fpkm_tracking)
			self.setOutput('isoforms_fpkm_tracking',isoforms_fpkm_tracking)
			self.setOutput('skipped_gtf',skipped_gtf)
			self.setOutput('transcripts_gtf',transcripts_gtf)
		else:
			self.setOutput('genes_fpkm_tracking',None)
			self.setOutput('isoforms_fpkm_tracking',None)
			self.setOutput('skipped_gtf',None)
			self.setOutput('transcripts_gtf',None)

	def call(self, *args):

		bamUpstream = args[0]

		self.setParamIO('bamInput',bamUpstream.getOutput('bamOutput'))

	def _singleRun(self,i):
		bamInput = self.getInputList('bamInput')
		gtfInput = self.getParamIO('gtfInput')
		#fragBiasCorrectInput = self.getInputList('fragBiasCorrectInput')
		outputDir = self.getParamIO('outputDir')
		print(os.path.join(Configure.getTmpDir(), 'assemblies.txt'))

		cmdline = [
				'cufflinks',
				'-p',str(self.getParam('threads')),
				self.getBoolParamCmd('-u','ismultiReadCorrect'),
				self.getBoolParamCmd('-N','isupperQuartileForm'),
				self.getBoolParamCmd('--total-hits-norm','istotalHitsNorm'),
				'-m',str(self.getParam('fragLenMean')),
				'-s',str(self.getParam('fragLenStdDev')),
				'-G',gtfInput,
				'-o',os.path.join(outputDir,'cufflinks_'+str(i)),
				bamInput[i],
				';',
				# 'echo', '"'+self.convertToRealPath(os.path.join(outputDir,'cufflinks_'+str(i),'transcripts.gtf')).split('.tmp')[1]+'" >>',
				'echo', '"'+os.path.join(outputDir,'cufflinks_'+str(i),'transcripts.gtf')+'" >>',
				self.getOutput("assembliesOutput")
				]

		result = self.callCmdline('V1', cmdline)
		f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
		f.write(result.stdout)
		f.close()

	def getMarkdownEN(self,):
		mdtext = """
### cufflinks Result
The cufflinks result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
Total map reads means that total number of reads mapped to genome
        """.format(mapRs = self.getOutput('stdOutput'))
		return mdtext