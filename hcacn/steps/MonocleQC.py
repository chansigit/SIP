# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from ..core import Step,Configure
import subprocess
import os

class MonocleQC(Step):
    """
    MonocleQC is a Step to filter genes and cells of single-cell dataset.
    The input file should be UMI counts data or mRNA counts data.
    See _init_ to initialize this step.
    > MonocleQC():_init_paramters
        matrixdata: str
    	A file(a table with txt format) path. The first row should be cells label, 
    	and the first should be genes name. Each row indicates a cell's data.
    outputpath: str
        A str indicates the name of appointed folder that saves outputs.You should
        build that folder in advance. The absolute path is also legel. 
    min_expression: int
        The expression threshold, default is 0.1. See more detailed information from 
        monocle package:"Sets the global expression detection threshold  to be used 
        with this CellDataSet. Counts how many cells each feature in a CellDataSet 
        object that are detectably expressed above a minimum threshold. Also counts
        the number of genes above this threshold are detectable in each cell."
    num_cells_expressed_threshold: int
        The thredshold of how many cell expressing a specific gene. Default is 10.
        The gene will be filtered oout if it expressed in less cells than this threshold.
    TotalmRNAs: float
        The threadshold of total_mRNAs expression of a individual cell. Default is 1e6.
        If the total_mRNA expression in the cell is greater than this threshold, the 
        cell will be filtered out.
    mean_expression_threshold: float
        The threshold of genes mean expression. Default is 0.1.If the gene's mean 
        expression is lower than this threshold, the gene will be filtered out.
    cmdParam: str or list of string
        current unsupported
    > MonocleQC()():_call_parameters
        Avaliabel upstream objects combinations:
            (DigitalExpression)
    """
    def __init__(self,
                 matrixdata = None,
                 outputpath = None,
                 min_expression = 0.1,
                 num_cells_expressed_threshold = 10,
                 TotalmRNAs = 1e6, 
                 mean_expression_threshold=0.1,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        """
        called by 'MonocleQC()'
        __init__(): Initialize the class with inputs, outputs and other parameters.
        Setting all parameter is the main target of this function.
        > Parameters
        matrixdata: str
            A file(a table with txt format) path. 
        outputpath: str
            The output directory of output results.
        min_expression: int
            The expression threshold, default is 0.1.
        num_cells_expressed_threshold: int
            The thredshold of how many cell expressing a specific gene. Default is 10.
        TotalmRNAs: float
            The threadshold of total_mRNAs expression of a individual cell. Default is 1e6.
        mean_expression_threshold: float
            The threshold of genes mean expression. Default is 0.1.
        cmdParam: str or list of string
            current unsupported
        """
       
        # set all input and output parameters
        self.setParamIO('matrixdata',matrixdata)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters
        self.setParam('min_expression',min_expression)
        self.setParam('num_cells_expressed_threshold',num_cells_expressed_threshold)
        self.setParam('TotalmRNAs',TotalmRNAs)
        self.setParam('mean_expression_threshold',mean_expression_threshold)
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
   
        self.setOutputDir1To1('density_Total_mRNAs',outputpath, 'density_Total_mRNAs','jpg','matrixdata')
        self.setOutputDir1To1('meanexpression_disersionemprical',outputpath,'meanexpression_disersionemprical','jpg','matrixdata')
        self.setOutputDir1To1('PCvariance', outputpath,'PCvariance','jpg','matrixdata')
        self.setOutputDir1To1('MonocleQCimage', outputpath,'MonocleQCimage','Rdata','matrixdata')
        self.setInputRscript('Rscript','MonocleQC.R')

    def call(self, *args):
        """
        called by Seurat()(object)
        """
        # the first object
        dropseqUpstream = args[0]      
        
        # set all required input parameters from upstream object
        self.setParamIO('matrixdata',dropseqUpstream.getOutput('dgeOutput'))
        
    def _multiRun(self,):
        matrixdata = self.getInput('matrixdata')
        min_expression = self.getParam('min_expression')
        num_cells_expressed_threshold = self.getParam('num_cells_expressed_threshold')
        TotalmRNAs = self.getParam('TotalmRNAs')
        mean_expression_threshold = self.getParam('mean_expression_threshold')
        Rscript = self.getInput('Rscript')
        cmdline = ['Rscript',
                    Rscript,
        			matrixdata,
                    str(min_expression),
                    str(num_cells_expressed_threshold),
                    str(TotalmRNAs),
                    str(mean_expression_threshold),
                    self.getParamIO('outputpath')
        			]
        print(''.join(cmdline))
        self.callCmdline('V1',cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## Monocle QC Result
```{{r, echo=FALSE, eval=FALSE}}
#don't run
MonocleQC(matrixdata='', outputpath='') 
```
### The Total_mRNAs~density Curve of all cells is shown below:

![]({density_Total_mRNAs})

### The picture below shows how variability (dispersion) in a gene's expression depends on the average expression across cells:
![]({meanexpression_disersionemprical})

### The picture below shows the variance explained by each component:
![]({PCvariance})

""".format(density_Total_mRNAs = self.getOutput('density_Total_mRNAs')[0],
                   meanexpression_disersionemprical = self.getOutput('meanexpression_disersionemprical')[0],
                   PCvariance = self.getOutput('PCvariance')[0])
        return mdtext
