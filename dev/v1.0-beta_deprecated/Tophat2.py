# coding: utf-8

from ..core import Step,Configure,Schedule
from .FastqDump import FastqDump
import os
class Tophat2(Step):
    def  __init__(self,
                  fastqInput1 = None,
                  fastqInput2 = None,
                  bt2Idx = None,
                  gtfInput = None,
                  outputDir = None,
                  threads = None,
                  cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        Configure.setRefSuffix('bt2Idx','.bt2.index/bt2_index',check=False)
        Configure.setRefSuffix('gtfInput','.gtf',check=False)
        # set all input and output parameters
        self.setParamIO('fastqInput1',fastqInput1)
        self.setParamIO('fastqInput2',fastqInput2)
        self.setParamIO('bt2Idx',bt2Idx)
        self.setParamIO('gtfInput',gtfInput)
        self.setParamIO('outputDir',outputDir)

        # call self.initIO()
        self.initIO()

        #set other parameters
        if threads is None:
            self.setParam('threads',Configure.getThreads())
        else:
            self.setParam('threads',threads)

    def impInitIO(self,):
        fastqInput1 = self.getParamIO('fastqInput1')
        fastqInput2 = self.getParamIO('fastqInput2')
        outputDir = self.getParamIO('outputDir')
        bt2Idx = self.getParamIO('bt2Idx')
        gtfInput = self.getParamIO('gtfInput')
        if outputDir is None:
            self.setParamIO('outputDir',Configure.getTmpDir())
        self.setOutput('stdOutput', os.path.join(Configure.getTmpDir(),'stdout.txt'))

        self.setInputDirOrFile('fastqInput1',fastqInput1)
        self.setInputDirOrFile('fastqInput2',fastqInput2)

        if bt2Idx is None:
            self.setParamIO('bt2Idx', Configure.getConfig('bt2Idx'))
            bt2Idx = self.getParamIO('bt2Idx')
        suffix = ['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']
        bt2IdxFiles = [ bt2Idx + s for s in suffix]
        self.setInput('bt2IdxFiles', bt2IdxFiles)

        if gtfInput is None:
            gtfInput=Configure.getConfig('gtfInput')
            self.setParamIO('gtfInput',gtfInput)
        self.setInput('gtfInput',gtfInput)

        if fastqInput1 is not None:
            self._setInputSize(len(self.getInputList('fastqInput1')))
            unmapped=list()
            accepted_hits=list()
            junctions=list()
            insertions=list()
            deletions=list()
            align_summary=list()
            prep_reads=list()
            for i in range(len(self.getInputList('fastqInput1'))):
                unmapped.append(os.path.join(outputDir, 'tophat_'+str(i),'unmapped.bam'))
                accepted_hits.append(os.path.join(outputDir, 'tophat_'+str(i),'accepted_hits.bam'))
                junctions.append(os.path.join(outputDir, 'tophat_'+str(i),'junctions.bed'))
                insertions.append(os.path.join(outputDir, 'tophat_'+str(i),'insertions.bed'))
                deletions.append(os.path.join(outputDir, 'tophat_'+str(i),'deletions.bed'))
                align_summary.append(os.path.join(outputDir, 'tophat_'+str(i),'align_summary.txt'))
                prep_reads.append(os.path.join(outputDir, 'tophat_'+str(i),'prep_reads.info'))
            self.setOutput('unmapped',unmapped)
            self.setOutput('samOutput',accepted_hits)###
            self.setOutput('junctions',junctions)
            self.setOutput('insertions',insertions)
            self.setOutput('deletions',deletions)
            self.setOutput('align_summary',align_summary)
            self.setOutput('prep_reads',prep_reads)
        else:
            self.setOutput('unmapped',None)
            self.setOutput('samOutput',None)###
            self.setOutput('junctions',None)
            self.setOutput('insertions',None)
            self.setOutput('deletions',None)
            self.setOutput('align_summary',None)
            self.setOutput('prep_reads',None)


    def call(self,*args):

        Upstream = args[0]

        if isinstance(Upstream,FastqDump):
            self.setParamIO('fastqInput1',Upstream.getOutput('fastqOutput1'))
            self.setParamIO('fastqInput2',Upstream.getOutput('fastqOutput2'))

    def _singleRun(self,i):
        gtfInput = self.getParamIO('gtfInput')
        fastqInput1 = self.getInputList('fastqInput1')
        fastqInput2 = self.getInputList('fastqInput2')
        outputDir = self.getParamIO('outputDir')
        bt2IdxFile = self.getParamIO('bt2Idx')

        cmdline = ['tophat2',
                  '-p',str(self.getParam('threads')),
                  '-G',gtfInput,
                  '-o',os.path.join(outputDir,'tophat_'+str(i)),
                  bt2IdxFile,
                  fastqInput1[i],
                  fastqInput2[i]]

        result = self.callCmdline('V2', cmdline)
        f = open(self.convertToRealPath(self.getOutput('stdOutput')),'ab+')
        f.write(result.stdout)
        f.close()
            
    def getMarkdownEN(self,):
        mdtext = """
### tophat2 Result
The tophat2 result is shown below:
```{{r, echo=FALSE}}
con <- file("{mapRs}", "r", blocking = FALSE)
readLines(con)
```
        """.format(mapRs = self.getOutput('stdOutput'))

        return mdtext