from ..core import Step,Configure,Schedule
import os

class Stringtie(Step):
	def __init__(self,
		bamInput = None,
		gtfInput = None,
		outputDir = None,
		threads = None,
		cmdParam = None,
		**kwargs
		):
		super(Step, self).__init__(cmdParam,**kwargs)
		Configure.setRefSuffix('gtfInput','.gtf',check=False)
		self.setParamIO('bamInput',bamInput)
		self.setParamIO('gtfInput',gtfInput)
		self.setParamIO('outputDir',outputDir)
		self.initIO()

		if threads is None:
			threads = Configure.getThreads()
		self.setParam('threads',threads)

	def impInitIO(self,):
		bamInput = self.getParamIO('bamInput')
		gtfInput = self.getParamIO('gtfInput')
		outputDir = self.getParamIO('outputDir')
		if outputDir is None:
			self.setParamIO('outputDir',Configure.getTmpDir())

		if gtfInput is None:
			gtfInput=Configure.getConfig('gtfInput')
			self.setParamIO('gtfInput',gtfInput)
		self.setInput('gtfInput',gtfInput)

		self.setInputDirOrFile('bamInput',bamInput)

		self.setOutput('assembliesOutput',os.path.join(Configure.getTmpDir(), 'assemblies.txt'))
		self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))
		self.setOutputDir1To1('transcripts_gtfOutput', outputDir, None, 'gtf', 'bamInput')

		if bamInput is not None:
			self._setInputSize(len(self.getInputList('bamInput')))

	def call(self, *args):

		bamUpstream = args[0]

		self.setParamIO('bamInput',bamUpstream.getOutput('bamOutput'))

	def _singleRun(self,i):
		bamInput = self.getInputList('bamInput')
		gtfInput = self.getParamIO('gtfInput')
		outputDir = self.getOutputList('transcripts_gtfOutput')

		cmdline = [
				'stringtie',
				bamInput[i],
				'-p',str(self.getParam('threads')),
				'-G',gtfInput,
				'-o',outputDir[i],
				';',
				# 'echo', '"'+self.convertToRealPath(os.path.join(outputDir,'cufflinks_'+str(i),'transcripts.gtf')).split('.tmp')[1]+'" >>',
				'echo', '"'+outputDir[i]+'" >>',
				self.getOutput("assembliesOutput")
				]

		result = self.callCmdline('V1', cmdline)
		f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
		f.write(result.stdout)
		f.close()

	def getMarkdownEN(self,):
		mdtext = """
### stringtie Result
The stringtie result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
        """.format(mapRs = self.getOutput('stdOutput'))
		return mdtext
