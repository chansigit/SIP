# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from ..core import Step,Configure
import os

class Monocle2QC(Step):
    """
    MonocleQC is a Step to filter genes and cells of single-cell dataset.
    The input file should be UMI counts data or mRNA counts data.
    See _init_ to initialize this step.
    > Monocle2QC():_init_paramters
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
                 featuredata = None,
                 phenodata = None,
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
        if featuredata:
            self.setParam('featuredata_flag',True)
            self.setParamIO('featuredata',featuredata)
        if phenodata:
            self.setParam('phenodata_flag',True)
            self.setParamIO('phenodata',phenodata)           
        # call self.initIO()
        self.initIO()
        #set other parameters
        self.setParam('min_expression',min_expression)
        self.setParam('num_cells_expressed_threshold',num_cells_expressed_threshold)
        self.setParam('TotalmRNAs',TotalmRNAs)
        self.setParam('mean_expression_threshold',mean_expression_threshold)
        #self._setMultiRun() ??
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
        matrixdata = self.getParamIO('matrixdata')        
        outputpath = self.getParamIO('outputpath')  
        featuredata_flag = self.getParam('featuredata_flag')
        phenodata_flag = self.getParam('phenodata_flag')

        if featuredata_flag:
            featuredata = self.getParamIO('featuredata')
            self.setInputDirOrFile('featuredata',featuredata)
            #
        if phenodata_flag:
            phenodata = self.getParamIO('phenodata')
            self.setInputDirOrFile('phenodata',phenodata)

        self.setInputDirOrFile('matrixdata',matrixdata)
        
        self.setInputRscript('Rscript','Monocle2QC.R')

        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath')  

        #MD5 
        def func1(basename):
            return basename+'/meanexpression_disersionemprical.jpg'
        def func2(basename):
            return basename+'/density_Total_mRNAs.jpg'     
        def func3(basename):
            return basename+'/PCvariance.jpg'
        def func4(basename):
            return basename+'/'+basename+'_cds.Rdata'
        
        self.setOutputDir1To1ByFunc('Monocle2QC_meanexpression_disersionemprical.jpg',outputpath, func1,'matrixdata')
        self.setOutputDir1To1ByFunc('Monocle2QC_density_Total_mRNAs.jpg',outputpath, func2,'matrixdata')
        self.setOutputDir1To1ByFunc('Monocle2QC_PCvariance.jpg',outputpath, func3,'matrixdata')
        self.setOutputDir1To1ByFunc('cds_file',outputpath, func4,'matrixdata')

        
        # set how many sample are there
        if matrixdata is not None:
            self._setInputSize(len(self.getInputList('matrixdata')))

    def call(self, *args):
        """
        called by Seurat()(object)
        """
        # the first object
        dropseqUpstream = args[0]      
        
        # set all required input parameters from upstream object
        #self.setParamIO('matrixdata',dropseqUpstream.getOutput('dgeOutput'))
        


    def _singleRun(self,i):
        # obtain all input and output dir list
        matrixdatas = self.getInputList('matrixdata')
        featuredata_flag = self.getParam('featuredata_flag')
        phenodata_flag = self.getParam('phenodata_flag')
        if featuredata_flag:
            featuredatas = self.getInputList('featuredata')
        else:
            featuredatas = ["None" for matrix in matrixdatas]
        if phenodata_flag:
            phenodatas = self.getInputList('phenodata')
        else:
            phenodatas = ["None" for matrix in matrixdatas]
        min_expression = self.getParam('min_expression')
        num_cells_expressed_threshold = self.getParam('num_cells_expressed_threshold')
        TotalmRNAs = self.getParam('TotalmRNAs')
        mean_expression_threshold = self.getParam('mean_expression_threshold')
        Rscript = self.getInput('Rscript')

        cds_files = self.getOutput('cds_file')
        cmdline = ['Rscript',
                    Rscript,
                    matrixdatas[i],
                    featuredatas[i],
                    phenodatas[i],
                    str(min_expression),
                    str(num_cells_expressed_threshold),
                    str(TotalmRNAs),
                    str(mean_expression_threshold),
                    os.path.dirname(cds_files[i])
                    ]
        self.callCmdline('V1', cmdline)


    def getMarkdownEN(self,):

        meanexpression_disersionemprical = self.getOutput('Monocle2QC_meanexpression_disersionemprical.jpg')
        density_Total_mRNAs = self.getOutput('Monocle2QC_density_Total_mRNAs.jpg')
        PCvariance = self.getOutput('Monocle2QC_PCvariance.jpg')

        meanexpression_disersionemprical_sen = ['***For %s***\n![meanexpression_disersionemprical Cluster](%s)\n'%(item.split("/")[-2],item) for item in meanexpression_disersionemprical]
        density_Total_mRNAs_sen = ['***For %s\n***\n![density_Total_mRNAs](%s)\n'%(item.split("/")[-2],item) for item in density_Total_mRNAs]
        PCvariance_sen = ['***For %s\n***\n![PCvariance](%s)\n'%(item.split("/")[-2],item)for item in PCvariance]
        
        meanexpression_disersionemprical_sen = "\n".join(meanexpression_disersionemprical_sen)
        density_Total_mRNAs_sen = "\n".join(density_Total_mRNAs_sen)
        PCvariance_sen = "\n".join(PCvariance_sen)


        mdtext = """
## Monocle2QC Usage

Monocle2QC('/path/to/cds.RData','/path/to/output_dir')  
    > Monocle2QC():_init_paramters
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
## Monocle2 Quality Control Result  
The Monocle2 Quality Control result is shown below:  
### meanexpression_disersionemprical jpg  
{meanexpression_disersionemprical}   

### density_Total_mRNAs jpg  

{density_Total_mRNAs} 
### PCvariance jpg  

{PCvariance} 



""".format(meanexpression_disersionemprical=meanexpression_disersionemprical_sen ,
    density_Total_mRNAs=density_Total_mRNAs_sen ,
    PCvariance =PCvariance_sen,
        )

    
        return mdtext
