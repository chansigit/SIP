
from ..core import Step,Configure
import subprocess
import os

class Cuffmerge(Step):
    Configure.setRefSuffix('faInput','.fa',check=False)
    Configure.setRefSuffix('gtfInput','.gtf',check=False)
    def __init__(self,
                 faInput = None,  
                 gtfInput = None,  
                 assembliesInput = None,
                 threads = None,
                 gtfOutputDir = None, 
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        
        self.setParamIO('faInput',faInput)
        self.setParamIO('gtfInput',gtfInput)
        self.setParamIO('assembliesInput',assembliesInput)
        self.setParamIO('gtfOutputDir',gtfOutputDir)

       
        self.initIO()
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads',threads)
            
        
    def impInitIO(self,):        
        faInput = self.getParamIO('faInput')
        gtfInput = self.getParamIO('gtfInput')
        assembliesInput = self.getParamIO('assembliesInput')
        gtfOutputDir = self.getParamIO('gtfOutputDir')
        if gtfOutputDir is None:
            self.setParamIO('gtfOutputDir',Configure.getTmpDir())
        self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))
        #set all input files        
        self.setInputDirOrFile('assembliesInput',assembliesInput) 
        
        if faInput is None:
            faInput=Configure.getConfig('faInput')
            self.setParamIO('faInput',faInput)
        self.setInput('faInput',faInput)

        if gtfInput is None:
            gtfInput=Configure.getConfig('gtfInput')
            self.setParamIO('gtfInput',gtfInput)
        self.setInput('gtfInput',gtfInput)

        if assembliesInput is not None:
            self._setInputSize(len(self.getInputList('assembliesInput')))
            merged_gtf=list()
            for i in range(len(self.getInputList('assembliesInput'))):
                merged_gtf.append(os.path.join(gtfOutputDir, 'cuffmerge_'+str(i),'merged.gtf'))
            self.setOutput('merged_gtf',merged_gtf)
        else:
            self.setOutput('merged_gtf',None)
        
    def call(self, *args):
        htseqUpstream = args[0]              
        self.setParamIO('assembliesInput',htseqUpstream.getOutput('assembliesOutput'))
        # self.setParamIO('gtfInput',htseqUpstream.getOutput('gtfOutput1'))
            
    def _singleRun(self, i):
        faInput = self.getParamIO('faInput')
        gtfInput = self.getParamIO('gtfInput')
        assembliesInput = self.getInputList('assembliesInput')
        gtfOutputDir = self.getParamIO('gtfOutputDir')
        cmdline = ['cuffmerge',
                    '-g', gtfInput,
                    '-s', faInput,
                    '-o', os.path.join(gtfOutputDir, 'cuffmerge_'+str(i)),
                    '-p', str(self.getParam('threads')),
                    assembliesInput[i]
                    ]
                    
        result = self.callCmdline('V2', cmdline)
        f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
        f.write(result.stdout)
        f.close()
            
    def getMarkdownEN(self,):
        mdtext = """
### cuffmerge Result
The cuffmerge result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
Total map reads means that total number of reads mapped to genome
        """.format(mapRs = self.getOutput('stdOutput'))

        return mdtext