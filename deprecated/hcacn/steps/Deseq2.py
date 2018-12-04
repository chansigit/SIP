# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
2018/3/30
"""
from ..core import Step,Configure
import subprocess
import os

class Deseq2(Step):
    def __init__(self,
                 matrixdata = None,
                 annotation = None,
                 outputpath = None,
                 padj = 0.05,
                 lgFDl = -1,
                 lgFDu = 1,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
       
        # set all input and output parameters
        self.setParamIO('matrixdata',matrixdata)
        self.setParamIO('annotation',annotation)
        self.setParamIO('outputpath',outputpath)
        # call self.initIO()
        self.initIO() 
        #set other parameters
        self.setParam('padj',padj)
        self.setParam('lgFDl',lgFDl)
        self.setParam('lgFDu',lgFDu)
        self._setMultiRun()

    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
        matrixdata = self.getParamIO('matrixdata')
        annotation = self.getParamIO('annotation')        
        outputpath = self.getParamIO('outputpath')  

        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath')  
        self.setInput('matrixdata', matrixdata)
        self.setInput('annotation', annotation)
        # create output file paths and set
        self.setOutputDir1To1('Deseq2Result',outputpath, 'Deseq2Result','csv','matrixdata')
        self.setOutputDir1To1('Summary',outputpath, 'Summary','txt','matrixdata')
        self.setOutputDir1To1('MAplot',outputpath, 'MAplot','jpg','matrixdata')
        self.setOutputDir1To1('DEheatmap',outputpath, 'DEheatmap','png','matrixdata')
        #set R script
        self.setInputRscript('Rscript', 'Deseq2.R')

    def call(self, *args):
        # the first object
        #dropseqUpstream = args[0]      
        # set all required input parameters from upstream object
        #self.setParamIO('matrixdata',dropseqUpstream.getOutput('dgeOutput'))
        pass

    def _multiRun(self,):
        matrixdata = self.getInput('matrixdata')
        annotation = self.getInput('annotation')
        padj = self.getParam('padj')
        lgFDl = self.getParam('lgFDl')
        lgFDu = self.getParam('lgFDu')
        Rscript = self.getInput('Rscript')
        cmdline = ['Rscript',
                    Rscript,
        			matrixdata,
                    annotation,
                    self.getParamIO('outputpath'),
                    str(padj),
                    str(lgFDl),
                    str(lgFDu)]
        print(''.join(cmdline))
        self.callCmdline('V1',cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## Gene Differential Expression Results
### Deseq2 have be done, and results can be found in the specific output folder and the folder named "result".
### The MA plot:
![]({MAplot})

### The heatmap of top 50 differential genes with high mean expression:
![]({DEheatmap})

""".format(MAplot = self.getOutput('MAplot')[0],
           DEheatmap = self.getOutput('DEheatmap')[0])
        return mdtext