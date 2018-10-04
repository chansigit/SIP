
# coding: utf-8
from ..core import Step,Configure,Schedule
import os
from .FastqDump import FastqDump
class FastQC(Step):

    def  __init__(self,

                  fastqInput=None,
                  fastqcOutputDir = None,
                  threads = None,
                  cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)

        # set all input and output parameters
        self.setParamIO('fastqInput',fastqInput)
        self.setParamIO('fastqcOutputDir',fastqcOutputDir)

        # call self.initIO()
        self.initIO()

        #set other parameters
        #self.setParam('isNoDiscordant', isNoDiscordant)
        #self.setParam('fileFormat',fileFormat)
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads',threads)

        print (self.params)

    def impInitIO(self,):

        # obtain all input and output parameters
        fastqInput = self.getParamIO('fastqInput')
        fastqcOutputDir = self.getParamIO('fastqcOutputDir')

        #set all input files
        self.setInputDirOrFile('fastqInput',fastqInput)


        # create output file paths and set
        if fastqcOutputDir is None:
            self.setParamIO('fastqcOutputDir',Configure.getTmpDir())
            fastqcOutputDir = self.getParamIO('fastqcOutputDir')

        self.setOutputDir1To1('fastqcOutput_html', fastqcOutputDir,None,'fastqc.html','fastqInput',sep='_')
        self.setOutputDir1To1('fastqcOutput_zip', fastqcOutputDir,None,'fastqc.zip','fastqInput',sep='_')
        
        


        # set how many sample are there
        if fastqInput is not None:
            self._setInputSize(len(self.getInputList('fastqInput')))


    def call(self,*args):

        Upstream = args[0]
        if isinstance(Upstream,FastqDump):
            fastqInput = Upstream.getOutput('fastqOutput1')
            fastqInput.extend( Upstream.getOutput('fastqOutput2'))
            self.setParamIO('fastqInput', fastqInput)
        # set all required input parameters from upstream object
        #上游可能為 “fastqInput1”,“fastqInput2”,“fastqOnput1”
        #self.setParamIO('fastqInput',Upstream.getOutput('fastqOutput1'))

        print("Call the UpStream Node, Not implementation")

        #other things

    def getMarkdownEN(self,):
        fastqcHtml = self.getOutput('fastqcOutput_html')
        fastqName = [ os.path.basename(item)[:-12] for item in fastqcHtml ]
        
        
        fastqName_R = "c(\"" + "\",\"".join(fastqName) + "\")"
        fastqcHtml_R = "c(\"" + "\",\"".join(fastqcHtml) + "\")"
        mdtext = """
## FastQC Usage

FastQC('/path/to/input_fastq',[fastq/sra],'/path/to/output_dir')  
Attention!  
* In this release, only .fastq format file can be setting as input!  

## FastQC Quality Control Result  
The FastQC Quality Control is shown below:  
```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
library(knitr)
library(kableExtra)
fastq_name <- {fastq_name}
fastqc_html <- {fastqc_report}
fq <- cbind(fastq_name, fastqc_html)
colnames(fq) <- c("Fastq Name", "Fastq Report")
kable(fq, "html") %>% kable_styling() %>% scroll_box(width = "1100px", height = "500px")
```

        
""".format(fastq_name = fastqName_R,fastqc_report= fastqcHtml_R )
 
        return mdtext
            
            
    def _singleRun(self,i):
        # obtain all input and output dir list
        fastqInput = self.getInputList('fastqInput')
        fastqcOutputDir = self.getParamIO('fastqcOutputDir')

        cmdline =['fastqc',
                   '-t',str(self.getParam('threads')),
                   '-o',fastqcOutputDir,
                   fastqInput[i]
                   ]
        self.callCmdline('V1', cmdline)
