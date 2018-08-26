# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
2018/3/29
"""

from ..core import Step,Configure
import subprocess
import os

class matrixRw(Step):

    def __init__(self,
                 matrixdata = None,
                 outputpath = None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
       
        # set all input and output parameters
        self.setParamIO('matrixdata',matrixdata)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters
        self._setMultiRun()
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
        matrixdata = self.getParamIO('matrixdata')        
        outputpath = self.getParamIO('outputpath')  

        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath')  
        self.setInput('matrixdata', matrixdata)
        # create output file paths and set
        self.setOutputDir1To1('processedMatrix',outputpath, 'processedMatrix','txt','matrixdata')
        self.setInputRscript('Rscript','matrixRw.R')


    def call(self, *args):
        # the first object
        #dropseqUpstream = args[0]      
        # set all required input parameters from upstream object
        #self.setParamIO('matrixdata',dropseqUpstream.getOutput('dgeOutput'))
        pass

    def _multiRun(self,):
        matrixdata = self.getInput('matrixdata')
        Rscript = self.getInput('Rscript')
        cmdline = ['Rscript',
                    Rscript,
        			matrixdata,
                    self.getParamIO('outputpath')
        			]
        print(''.join(cmdline))
        self.callCmdline('V1',cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## Pre-process Results
### The pre-process have be done.
"""
        return mdtext