from ..core import Step,Configure
import os

class AdapterRemoval(Step):
    def __init__(self,
                 fastqInput1 = None, 
                 fastqInput2 = None,                  
                 fastqOutputDir1 = None, 
                 fastqOutputDir2 = None,
                 adapter1 = None,
                 adapter2 = None,
                 threads = None,
                 cmdParamFindAdapter = None,
                 cmdParam = None, 
                 **kwargs):
        super(Step, self).__init__(cmdParam = [cmdParamFindAdapter, cmdParam],**kwargs)
        """
        called by 'AdapterRemoval()'
        __init__(): Initialize the class with inputs, outputs and other parameters.
        Setting all parameter is the main target of this function.
        """
        # set all input and output parameters
        self.setParamIO('fastqInput1',fastqInput1)
        self.setParamIO('fastqInput2',fastqInput2)
        self.setParamIO('fastqOutputDir1',fastqOutputDir1)
        self.setParamIO('fastqOutputDir2',fastqOutputDir2)    
       
        # call self.initIO()
        self.initIO()
            
        #set other parameters
        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads',threads)
        if adapter1 is None:
            self.setParam('adapter1',None)
        elif len(adapter1) > 1:
            self.setParam('adapter1',adapter1)
        elif os.path.exists(adapter1):
            self.setParam('adapter1',self.getListInFile())
        else:
            self.setParam('adapter1',adapter1)
        
        if adapter2 is None:
            self.setParam('adapter2',None)
        elif len(adapter1) > 1:
            self.setParam('adapter2',adapter2)
        elif os.path.exists(adapter2):
            self.setParam('adapter2',self.getListInFile())
        else:
            self.setParam('adapter2',adapter2)
        
        self.adapter1 = {}
        self.adapter2 = {}

        
    def impInitIO(self,): 
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters
        fastqInput1 = self.getParamIO('fastqInput1')
        fastqInput2 = self.getParamIO('fastqInput2')
        fastqOutputDir1 = self.getParamIO('fastqOutputDir1')
        fastqOutputDir2 = self.getParamIO('fastqOutputDir2')


        
        #set all input files
        self.setInputDirOrFile('fastqInput1',fastqInput1)       
        self.setInputDirOrFile('fastqInput2',fastqInput2)

        
        # create output file paths and set

        if fastqOutputDir1 is None:
            self.setParamIO('fastqOutputDir1',Configure.getTmpDir())
        if fastqOutputDir2 is None:
            self.setParamIO('fastqOutputDir2',Configure.getTmpDir())
        self.setOutputDir1To1('fastqOutput1',fastqOutputDir1, None, 'fastq','fastqInput1')       
        self.setOutputDir1To1('fastqOutput2',fastqOutputDir2, None, 'fastq','fastqInput2')
        self.setOutputDir1To1('adapterOutput', None, None, 'adapter.txt', 'fastqInput1')
        self.setOutputDir1To1('settingsOutput', None, None,'settings', 'fastqInput1')
        
        
        # set how many sample are there
        if fastqInput1 is not None:
            self._setInputSize(len(self.getInputList('fastqInput1')))
        
    def call(self, *args):
        """
        called by AdapterRemoval()(upstreamObj)
        """
        # the first object
        fastqUpstream = args[0]

        self.setParamIO('fastqInput1', fastqUpstream.getOutput('fastqOutput1'))
        self.setParamIO('fastqInput2', fastqUpstream.getOutput('fastqOutput2'))
            
    def _singleRun(self, i):
        """
        create and execute the command line        
        i is the No. of the sample
        """
        #get all input and output
        fastqInput1 = self.getInputList('fastqInput1')
        fastqInput2 = self.getInputList('fastqInput2')
        fastqOutput1 = self.getOutputList('fastqOutput1')
        fastqOutput2 = self.getOutputList('fastqOutput2')
        adapterOutput = self.getOutputList('adapterOutput')
        settingsOutput = self.getOutputList('settingsOutput')
        
        adapter1 = self.getParam('adapter1')
        if adapter1 is not None:
            adapter1 = adapter1[i]
        adapter2 = self.getParam('adapter2')
        if adapter2 is not None:
            adapter2 = adapter2[i]
        
        if adapter1 is None or adapter2 is None:
            #combine the command line
            cmdline = [
                    'AdapterRemoval',
                    '--file1',fastqInput1[i],
                    '--file2',fastqInput2[i],
                    '--identify-adapters',
                    '--threads',str(self.getParam('threads')),
                    '| grep Consensus: >',
                    adapterOutput[i],
                    ]
            #run commandline
            #self.callCmdline('V1',cmdline,stdoutToLog = False)
            self.callCmdline('V1',cmdline)
            #print(self.convertToRealPath(adapterOutput[i]))
            result = self.getListInFile(self.convertToRealPath(adapterOutput[i]))
            adapter1 = self.adapter1[str(i)] = result[0].split()[1]
            adapter2 = self.adapter2[str(i)] = result[1].split()[1]
        #combine the command line   
        cmdline  = [#'/root/software/adapterremoval/build/AdapterRemoval',
                'AdapterRemoval',
                '--file1',fastqInput1[i],
                '--file2',fastqInput2[i],
                '--adapter1',adapter1,
                '--adapter2',adapter2,
                '--output1',fastqOutput1[i],
                '--output2',fastqOutput2[i],
                '--threads',str(self.getParam('threads')),                
                '--basename', settingsOutput[i][0:-9]]  
        #run commandline
        self.callCmdline('V1',cmdline)
        
    def getMarkdownEN(self,):
        settingsOutput = self.getOutput('settingsOutput')
        settingsOutput_path = []

        for st in settingsOutput:
            settingsOutput_path.append("\"" + st + "\"")
        settingsOutput_path_str = ",".join([i for i in settingsOutput_path])
        fastqc1_path_str = "c(" + settingsOutput_path_str + ")"


        mdtext = """
## Adapter Removal Result
```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
library(knitr)
library(kableExtra)
st <- {settings}
st <- as.data.frame(st)
colnames(st) <- c("Adapter Information File Path")
```

The AdapterRemoval report path for every sample are as follows:

```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
kable(st, "html") %>% kable_styling() %>% scroll_box(width = "800px", height = "400px")
```

        """.format(settings=fastqc1_path_str)
        return mdtext


