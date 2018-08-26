

from ..core import Step,Configure
from ..steps import SingleCellExperiment
import os
class SC3_DE(Step):
    def __init__(self,
                 sceInput = None,
                 outputpath = None,
                 cmdParam=None,
                 **kwargs):
        """
        SC3_DE is XXX. Need to use this 
        Step as the downstream of SingleCellExperiment.
        >SC3_DE():_init_parameters
            sce: str
            The R workspace saved by SingleCellExperiment Step.
            outputpath: str
            A str indicates the name of appointed folder that saves outputs.You should
            build that folder in advance. The absolute path is also legel. 
            cluster_num: int
            The number of clusters.
            set to 0 will auto estimate the cluster number
            cmdParam: str or list of string
            current unsupported
        >SC3_DE()():_call_parameters
            Avaliabel upstream objects combinations:
            (SingleCellExperiment)
        """
        super(Step, self).__init__(cmdParam,**kwargs)
        
        

        # set all input and output parameters
        self.setParamIO('sceInput',sceInput)
        self.setParamIO('outputpath',outputpath) 
        # call self.initIO()
        self.initIO()
        #set other parameters



        #self._setMultiRun()
        
    def impInitIO(self,):
        """
        This function is to initialize 
        all of the input and output files from the io parameters set in __init__() 
        """
        # obtain all input and output parameters        
               
        outputpath = self.getParamIO('outputpath')  
        #set the input file
        if outputpath is None:
            self.setParamIO('outputpath',Configure.getTmpDir()) 
            outputpath = self.getParamIO('outputpath') 

        sceInput = self.getParamIO('sceInput') 
        #set all input files
        self.setInputDirOrFile('sceInput',sceInput)

        # create output file paths and set
        #self.setOutputDir1To1ByFunc('sceOutput',outputpath,func,"matrix_file")
        #self.setOutputDir1To1('sc3OutputFolder',outputpath,None,"_folder","sceInput",sep='')
        def func1(basename):
            return basename+'/De_genes.table'
        def func2(basename):
            return basename+'/Markers_genes.table'          
        self.setOutputDir1To1ByFunc('sc3Output_De_genes.table',outputpath, func1,'sceInput')
        self.setOutputDir1To1ByFunc('sc3Output_Markers_genes.table',outputpath, func2,'sceInput')

        # Rscripts
        self.setInputRscript('Rscript','SC3.R')

        if sceInput is not None:
            self._setInputSize(len(self.getInputList('sceInput')))
    def call(self,*args):

        Upstream = args[0]
        if isinstance(Upstream,SingleCellExperiment):
            self.setParamIO('sceInput', Upstream.getOutput('sceOutput'))


    def getMarkdownEN(self,):

        Result_table = self.getOutput('sc3Output_De_genes.table')
        list_name = [  item.split("/")[-2] for item in Result_table]
        list_name = "c(\"" + "\",\"".join(list_name) + "\")"
        list_table = "c(\"" + "\",\"".join(Result_table) + "\")"
        # print(list_table)
        # print(list_name)
        mdtext = """
## SC3_DE Usage

SC3_DE('/path/to/sce.RData','/path/to/output_dir')  
To find differential genes,Annotation file to build SingleCellExpeiment Object is needed!! the column name of labels in Annotation file must be "cell_type"
## SC3 Differential Expression Result  
The SC3 Differential Expression result is shown below:  
### Differential Expression Genes Result  

```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
library(knitr)
library(kableExtra)
list_name <- {list_name}
list_table <- {list_table}
for (i in 1:length(list_name))
{{
    #show(i)
    de_genes <- read.table(list_table[[i]])
    name <- list_name[[i]]
    colnames(de_genes) <- paste(name,"Pvalue",sep=".")
    col <- rownames(de_genes)[1:10]
    value <- de_genes[1:10,]
    data <- cbind(col, value)
    colnames(data) <- c(name, "Pvalue")
    #kable(data, "html") %>% kable_styling() %>% scroll_box(width = "1100px", height = "500px")
    print(data)

}}

```



  
""".format(  list_name =list_name, list_table =list_table)


        return mdtext
            
            
    def _singleRun(self,i):
        # obtain all input and output dir list
        sceInputs = self.getInputList('sceInput')
        sc3OutputTables = self.getOutput('sc3Output_De_genes.table')
        os.path.dirname
        Rscript = self.getInput('Rscript')
        cmdline =['Rscript',
                  Rscript,
                   sceInputs[i],
                   str(0),   # need no anything about cluster numbers
                   os.path.dirname(sc3OutputTables[i]),
                   'de'
                   ]
        self.callCmdline('V1', cmdline)