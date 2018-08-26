# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from ..core import Step,Configure
import subprocess
import os

class MonocleDC(Step):
    def __init__(self,
                 imageRdata = None,
                 outputpath = None,
                 num_PCA = 10,
                 cluster_num = 6,
                 cmdParam=None,
                 **kwargs):
        """
        MonocleDC is a Step to apply dimention reduction by T-SNE 
        and cluster cells by density peak clustering algorithm. Recommend to use this 
        Step as the downstream of MonocleQC, or it may lead to some errors.
        >MonocleDC():_init_parameters
            imageRdata: str
            The R workspace saved by MonocleQC Step.
            outputpath: str
            A str indicates the name of appointed folder that saves outputs.You should
            build that folder in advance. The absolute path is also legel. 
            num_PCA: int
            The number of pinciple components used in T-SNE dimention reduction.Default is 10.
            cluster_num: int
            The number of clusters.Default is 6.
            cmdParam: str or list of string
            current unsupported
        >MonocleDC()():_call_parameters
            Avaliabel upstream objects combinations:
            (MonocleQC)
        """
        super(Step, self).__init__(cmdParam,**kwargs)
        """
        called by 'MonocleDC()'
        __init__(): Initialize the class with inputs, outputs and other parameters.
        Setting all parameter is the main target of this function.
        >Parameters
        imageRdata: str
            The R workspace saved by MonocleQC Step.
        outputpath: str
            A str indicates the name of appointed folder that saves outputs.
        num_PCA: int
            The number of pinciple components used in T-SNE dimention reduction.Default is 10.
            cluster_num: int
        The number of clusters.
        cmdParam: str or list of string
            current unsupported
        """
        
        #Configure.enableDocker(False)
        # set all input and output parameters
        self.setParamIO('imageRdata',imageRdata)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters

        self.setParam('num_PCA',num_PCA)
        self.setParam('cluster_num',cluster_num)
        self._setMultiRun()
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
        imageRdata = self.getParamIO('imageRdata')        
        outputpath = self.getParamIO('outputpath')  

        
        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath')  
        self.setInput('imageRdata', imageRdata)
        # create output file paths and set
        self.setOutputDir1To1('densitypeak_cluster',outputpath, 'densitypeak_cluster','jpg','imageRdata')
        self.setInputRscript('Rscript','MonocleDC.R')
       
    def call(self, *args):
        """
        called by Seurat()(object)
        """
        # the first object
        MonocleQCupstream = args[0]      
        
        # set all required input parameters from upstream object
        self.setParamIO('imageRdata',MonocleQCupstream.getOutput('MonocleQCimage'))
        #print(MonocleQCupstream.getOutput('MonocleQCimage'))
        
    def _multiRun(self,):
        imageRdata = self.getInput('imageRdata')
        num_PCA = self.getParam('num_PCA')
        cluster_num = self.getParam('cluster_num')
        Rscript = self.getInput('Rscript')
        cmdline = ['Rscript',
                    Rscript,
        			imageRdata[0],
                    str(num_PCA),
                    str(cluster_num),
                    self.getParamIO('outputpath')
        			]
        print(''.join(cmdline))
        self.callCmdline('V1',cmdline)
    def getMarkdownEN(self,):
        mdtext="""
## Monocle dimention reduction and clustering results
```{{r, echo=FALSE,  eval=FALSE}}
#don't run
MonocleDC(outputpath='outputresult')(MonocleQC_result)
```
### Select the number of principle components used in dimension reduction, monocle use T-SNE to project the high-dimentional points into a low-dimentional space. Density-peak cluster method is used in clustering.
![]({densitypeak_cluster})

        """.format(densitypeak_cluster = self.getOutput('densitypeak_cluster')[0])
        return mdtext
