# -*- coding: utf-8 -*-
"""
Created on Tus Mar  7 12:28:52 2018

@author: ShengquanChen
"""

from ..core import Step,Configure
import subprocess
import os

class HTSeq_sam2count(Step):
    Configure.setRefSuffix('gtfInput','.gtf',check=False)
    def __init__(self,
                 samInput = None,  
                 gtfInput = None,             
                 countOutputDir = None, 
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        self.setParamIO('samInput',samInput)
        self.setParamIO('gtfInput',gtfInput)
        self.setParamIO('countOutputDir',countOutputDir)

       
        self.initIO()
            
        
    def impInitIO(self,):        
        samInput = self.getParamIO('samInput')
        gtfInput = self.getParamIO('gtfInput')
        countOutputDir = self.getParamIO('countOutputDir')
        if countOutputDir is None:
            self.setParamIO('countOutputDir',Configure.getTmpDir())
        self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))

        self.setInputDirOrFile('samInput',samInput) 
        if gtfInput is None:
            gtfInput = Configure.getConfig('gtfInput')
            self.setParamIO('gtfInput',gtfInput)
        self.setInput('gtfInput',gtfInput)
       
        self.setOutputDir1To1('countOutput', countOutputDir, None, 'count','samInput') 
        
        if samInput is not None:
            self._setInputSize(len(self.getInputList('samInput')))
        
    def call(self, *args):
        samUpstream = args[0]      
        
        self.setParamIO('samInput', samUpstream.getOutput('samOutput'))
            
    def _singleRun(self, i):
        samInput = self.getInputList('samInput')
        gtfInput = self.getParamIO('gtfInput')
        countOutput = self.getOutputList('countOutput')
        cmdline = ['python -m HTSeq.scripts.count', '-s no',
                    samInput[i],
                    gtfInput,
                    '>',
                    countOutput[i]
                    ]
        result = self.callCmdline('V1', cmdline)
        f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
        f.write(result.stdout)
        f.close()
            
    def getMarkdownEN(self,):
        mdtext = """
### HTSeq Result
The HTSeq result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
        """.format(mapRs = self.getOutput('stdOutput'))
        return mdtext