# coding: utf-8

from ..core import Step,Configure,Schedule
from .FastqDump import FastqDump
import os
class Hisat2(Step):
    Configure.setRefSuffix('ht2Idx','.ht2.index/genome',check=False)
    def  __init__(self,
                  fastqInput1 = None,
                  fastqInput2 = None,
                  ht2Idx = None,
                  samOutputDir = None,
                  threads = None,
                  cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)

        # set all input and output parameters
        #self.setParamIO('fastqInput1',fastqInput1)
        self.setParamIO('fastqInput1',fastqInput1)
        self.setParamIO('fastqInput2',fastqInput2)
        self.setParamIO('ht2Idx',ht2Idx)
        self.setParamIO('samOutputDir',samOutputDir)

        # call self.initIO()
        self.initIO()

        #set other parameters
        #self.setParam('isNoDiscordant', isNoDiscordant)
        if threads is None:
            self.setParam('threads',Configure.getThreads())
        else:
            self.setParam('threads',threads)

    def impInitIO(self,):

        # obtain all input and output parameters
        #fastqInput1 = self.getParamIO('fastqInput1')
        fastqInput1 = self.getParamIO('fastqInput1')
        fastqInput2 = self.getParamIO('fastqInput2')
        samOutputDir = self.getParamIO('samOutputDir')
        ht2Idx = self.getParamIO('ht2Idx')
        if samOutputDir is None:
            self.setParamIO('samOutputDir',Configure.getTmpDir())
        self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))

        #set all input files
        self.setInputDirOrFile('fastqInput1',fastqInput1)
        self.setInputDirOrFile('fastqInput2',fastqInput2)

        if ht2Idx is None:
            self.setParamIO('ht2Idx', Configure.getConfig('ht2Idx'))
            ht2Idx = self.getParamIO('ht2Idx')
        suffix = ['.1.ht2','.2.ht2','.3.ht2','.4.ht2','.5.ht2','.6.ht2','.7.ht2','.8.ht2']
        ht2IdxFiles = [ ht2Idx + s for s in suffix ]
        self.setInput('ht2IdxFiles', ht2IdxFiles)

        # create output file paths and set
        self.setOutputDir1To1('samOutput',samOutputDir,None,'sam','fastqInput1')

        # set how many sample are there
        if fastqInput1 is not None:
            self._setInputSize(len(self.getInputList('fastqInput1')))


    def call(self,*args):

        Upstream = args[0]

        # set all required input parameters from upstream object
        if isinstance(Upstream,FastqDump):
            self.setParamIO('fastqInput1',Upstream.getOutput('fastqOutput1'))
            self.setParamIO('fastqInput2',Upstream.getOutput('fastqOutput2'))

    def _singleRun(self,i):
        # obtain all input and output dir list
        fastqInput1 = self.getInputList('fastqInput1')
        fastqInput2 = self.getInputList('fastqInput2')
        samOutput = self.getOutputList('samOutput')

        ht2IdxFile = self.getParamIO('ht2Idx')

        cmdline = ['/root/software/hisat2-2.1.0/hisat2',
                  '-p',str(self.getParam('threads')),
                  '-x',ht2IdxFile,
                  '-1', fastqInput1[i],
                  '-2',fastqInput2[i],
                  '-S',samOutput[i]]

        result = self.callCmdline('V1', cmdline)
        f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
        f.write(result.stdout)
        f.close()
            
    def getMarkdownEN(self,):
        mdtext = """
### hisat Result
The hisat result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
Total map reads means that total number of reads mapped to genome
        """.format(mapRs = self.getOutput('stdOutput'))

        return mdtext