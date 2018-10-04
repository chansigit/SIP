from ..core import Step, Configure
from .Cellranger import Cellranger
import os

class Seuratpreprocessing(Step):
    """
        Seuratpreprocessing is a Step to preprocess reads matrix to Seurat object.
        See __init__ to initialize this step.
        > Seuratpreprocessing(): __init__  parameters
            inputdir: str or str list
                a directory contain all of matrix files, for example: /home/data/matrix
            outputdir: str
                outputdir path of all the results, default: ~/step_01_Seuratpreprocessing/
            rscript: str
                rscript file path, for example: /home/data/Seuratpreprocessing.R
            datatype: str
                datatype of input files, for example: '10x'
            cmdParam: str or list of string
                current unsupported
         > Seuratpreprocessing()(): __call__ parameters
            Available upstream objects combinations:
            (Celllranger)
            (to be added)
        """
    def __init__(self,
                 outputdir=None,
                 # barcodes=None,
                 # genes=None,
                 inputdir = None,
                 rscript=None,
                 datatype = None,
                 # matrix=None,
                 # densematrix = None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)
        # if densematrix or sparsematrix?
        # if matrix!=None and barcodes!=None and genes!=None:
        # self.setParamIO('matrix', matrix)
        # self.setParamIO('barcodes', barcodes)
        # self.setParamIO('genes', genes)
        # else:
        #     self.setParamIO('densematrix', densematrix)
        self.setParamIO('inputdir', inputdir)
        self.setParamIO('outputdir', outputdir)
        self.setParam('datatype', datatype)

        self.initIO()
        self._setMultiRun()

    def impInitIO(self, ):
        # matrix = self.getParamIO('matrix')
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        outputdir = self.getParamIO('outputdir')
        inputdir = self.getParamIO('inputdir')
        datatype = self.getParam('datatype')

        self.setInputRscript('rscript', '../rscript/Seuratpreprocessing.R')
        # if outputdir is None, os will error
        if outputdir is None:
           self.setParamIO('outputdir', Configure.getTmpDir())
           outputdir = self.getParamIO('outputdir')

        # set output/input paths
        # if inputdir is None, os will error
        if inputdir is not None:
            if datatype is'10x':
                self.setInputDirOrFile('barcodes', os.path.join(inputdir, 'barcodes.tsv'))
                self.setInputDirOrFile('genes', os.path.join(inputdir, 'genes.tsv'))
                self.setInputDirOrFile('matrix', os.path.join(inputdir, 'matrix.mtx'))
            else:
                self.setInputDirOrFile('matrix', inputdir)

            self.setOutputDirNTo1('variableGenes', os.path.join(outputdir, 'variableGenes.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('Rdata', os.path.join(outputdir, 'Preprocessing.Rdata'), '', 'barcodes')
            self.setOutputDirNTo1('violinplot', os.path.join(outputdir, 'violinplot.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('geneplot', os.path.join(outputdir, 'geneplot.jpeg'), '', 'barcodes')
            # self.setOutputDirNTo1('variableGenes', os.path.join(outputdir, 'variableGenes.jpeg'), '', 'barcodes')
            # self.setOutputDirNTo1('Elbowplot', os.path.join(outputdir, 'Elbowplot.jpeg'), '', 'barcodes')
            # self.setOutputDirNTo1('TSNEplot', os.path.join(outputdir, 'TSNEplot.jpeg'), '', 'barcodes')

    def call(self, *args):
        # the first object
        Upstream = args[0]
        # set all required input parameters from upstream objects
        if isinstance(Upstream, Cellranger):
            datatype = '10x'
            self.setParam('datatype', datatype)
        else:
            self.setParam('datatype', 'densemtx')
        self.setParamIO('inputdir', Upstream.getParamIO('finaldir'))
        # self.setParamIO('genes', cellrangerUpstream.getOutput('genes'))
        # self.setParamIO('matrix', cellrangerUpstream.getOutput('matrix'))

        # print(self.cmdline)

    def _multiRun(self, ):
        # get input parameters
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        # matrix = self.getParamIO('matrix')
        inputdir = self.getParamIO('inputdir')
        rscript = self.getInput('rscript')
        datatype = self.getParam('datatype')
        # get output parameters
        outputdir = self.getParamIO('outputdir')

        #
        # replace file path
        # inputdir = inputdir.replace('/home/cfeng', '/data')
        # # matrix = matrix.replace('/home/cfeng', '/data')
        # rscript = rscript.replace('/home/cfeng', '/data')
        # # barcodes = barcodes.replace('/home/cfeng', '/data')
        #
        # outputdir = outputdir.replace('/home/cfeng', '/data')
        cmdline = ['Rscript',
                   rscript,
                   inputdir,
                   datatype,
                   outputdir]
        print(cmdline)
        self.callCmdline(cmdline=cmdline, dockerVersion='V1')

    def getMarkdownEN(self,):
        rmd = """

## Data Preprocessing
`Preobject = Seuratpreprocessing(inputdir, outputdir, datatype=['10x','densemtx'], rscript)(upstream)`

You can set threshholds for downstream analysis based on below two pictures,
here we visualize gene and molecule counts, plot their relationship, thus you can exclude cells with a clear outlier number of genes detected as potential multiplets.

### Violinplot:
Violin plot is a method of plotting numeric data. It is similar to box plot with a rotated kernel density plot on each side.
![]({vioimagepath})

### Geneplot:
Geneplot here shows the relations between nGenes and nUMIs to determine the suitable threshold.
![]({geneimagepath})

### Select variable genes
Seurat calculates highly variable genes and focuses on these for downstream analysis.
![]({variableimagepath})
        """.format(vioimagepath=self.getOutput('violinplot'),
                   geneimagepath = self.getOutput('geneplot'),
                   variableimagepath=self.getOutput('variableGenes'))
        return rmd

