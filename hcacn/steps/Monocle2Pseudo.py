# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
"""

from ..core import Step,Configure
from .Monocle2QC import Monocle2QC
import subprocess
import os

class Monocle2Pseudo(Step):
    def __init__(self,
                 cds_file = None,
                 outputpath = None,
                 # num_PCA = 10,
                 # cluster_num = 6,
                 cmdParam=None,
                 **kwargs):
        """
        MonocleDC is a Step to apply dimention reduction by T-SNE 
        and cluster cells by density peak clustering algorithm. Recommend to use this 
        Step as the downstream of MonocleQC, or it may lead to some errors.
        >MonocleDC():_init_parameters
            cds_file: str
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
        cds_file: str
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
        self.setParamIO('cds_file',cds_file)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters

        # self.setParam('num_PCA',num_PCA)
        # self.setParam('cluster_num',cluster_num)

        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
        cds_file = self.getParamIO('cds_file')        
        outputpath = self.getParamIO('outputpath')  

        
        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath')  

        self.setInputDirOrFile('cds_file',cds_file)
        
        self.setInputRscript('Rscript','Monocle2Pseudo.R')           

        def func1(basename):
            return basename+'/ordering_genes.jpg'
        def func2(basename):
            return basename+'/cell_clusters.jpg'     
        def func3(basename):
            return basename+'/cell_trajectory_by_groups.jpg'
        def func4(basename):
            return basename+'/cell_trajectory_by_Pseudotime.jpg'
        def func5(basename):
            return basename+'/cell_trajectory_by_Cluster.jpg'     
        def func6(basename):
            return basename+'/cell_trajectory_by_State.jpg'
        def func7(basename):
            return basename+'/cds_Pseudo.Rdata'
        
        self.setOutputDir1To1ByFunc('Monocle2Pseudo_ordering_genes.jpg',outputpath, func1,'cds_file')
        self.setOutputDir1To1ByFunc('Monocle2Pseudo_cell_clusters.jpg',outputpath, func2,'cds_file')
        self.setOutputDir1To1ByFunc('Monocle2Pseudo_cell_trajectory_by_groups.jpg',outputpath, func3,'cds_file')
        self.setOutputDir1To1ByFunc('Monocle2Pseudo_cell_trajectory_by_Pseudotime.jpg',outputpath, func4,'cds_file')
        self.setOutputDir1To1ByFunc('Monocle2Pseudo_cell_trajectory_by_Cluster.jpg',outputpath, func5,'cds_file')
        self.setOutputDir1To1ByFunc('Monocle2Pseudo_cell_trajectory_by_State.jpg',outputpath, func6,'cds_file')
        self.setOutputDir1To1ByFunc('cds_Pseudo_file',outputpath, func7,'cds_file')
        # set how many sample are there
        if cds_file is not None:
            self._setInputSize(len(self.getInputList('cds_file')))
    def call(self, *args):
        """
        called by Seurat()(object)
        """
        # the first object
        MonocleQCupstream = args[0]      
        
        if isinstance(MonocleQCupstream,Monocle2QC):
            # set all required input parameters from upstream object
            self.setParamIO('cds_file',MonocleQCupstream.getOutput('cds_file'))
        #print(MonocleQCupstream.getOutput('MonocleQCimage'))
 
    def _singleRun(self,i):
        # obtain all input and output dir list
        cds_files = self.getInputList('cds_file')
        
        Rscript = self.getInput('Rscript')
        cds_Pseudo_files = self.getOutput('cds_Pseudo_file')
        print(cds_Pseudo_files)
        cmdline = ['Rscript',
                    Rscript,
                    cds_files[i],
                    os.path.dirname(cds_Pseudo_files[i])
                    ]
        self.callCmdline('V1', cmdline)       

    def getMarkdownEN(self,):
        cell_trajectory_by_groups = self.getOutput('Monocle2Pseudo_cell_trajectory_by_groups.jpg')
        cell_trajectory_by_Pseudotime = self.getOutput('Monocle2Pseudo_cell_trajectory_by_Pseudotime.jpg')
        cell_trajectory_by_Cluster = self.getOutput('Monocle2Pseudo_cell_trajectory_by_Cluster.jpg')
        cell_trajectory_by_State = self.getOutput('Monocle2Pseudo_cell_trajectory_by_State.jpg')

        cell_trajectory_by_groups_sen = ['***For %s***\n![cell_trajectory_by_groups Cluster](%s)\n'%(item.split("/")[-2],item) for item in cell_trajectory_by_groups]
        cell_trajectory_by_Pseudotime_sen = ['***For %s\n***\n![cell_trajectory_by_Pseudotime](%s)\n'%(item.split("/")[-2],item) for item in cell_trajectory_by_Pseudotime]
        cell_trajectory_by_Cluster_sen = ['***For %s\n***\n![cell_trajectory_by_Cluster](%s)\n'%(item.split("/")[-2],item)for item in cell_trajectory_by_Cluster]
        cell_trajectory_by_State_sen = ['***For %s\n***\n![cell_trajectory_by_State](%s)\n'%(item.split("/")[-2],item)for item in cell_trajectory_by_State]
        
        cell_trajectory_by_groups_sen = "\n".join(cell_trajectory_by_groups_sen)
        cell_trajectory_by_Pseudotime_sen = "\n".join(cell_trajectory_by_Pseudotime_sen)
        cell_trajectory_by_Cluster_sen = "\n".join(cell_trajectory_by_Cluster_sen)
        cell_trajectory_by_State_sen = "\n".join(cell_trajectory_by_State_sen)


        mdtext = """
## Monocle2Pseudo Usage

Monocle2Pseudo('/path/to/cds.RData','/path/to/output_dir')  
    
## Monocle2 Pseudotime Result  
The Monocle2 Pseudotime result is shown below:  
### cell_trajectory_by_groups jpg  
{cell_trajectory_by_groups}   

### cell_trajectory_by_Pseudotime jpg  

{cell_trajectory_by_Pseudotime} 
### cell_trajectory_by_Cluster jpg  

{cell_trajectory_by_Cluster} 

### cell_trajectory_by_State jpg  

{cell_trajectory_by_State} 

""".format(cell_trajectory_by_groups=cell_trajectory_by_groups_sen ,
    cell_trajectory_by_Pseudotime=cell_trajectory_by_Pseudotime_sen ,
    cell_trajectory_by_Cluster =cell_trajectory_by_Cluster_sen,
    cell_trajectory_by_State =cell_trajectory_by_State_sen,
        )


        return mdtext
