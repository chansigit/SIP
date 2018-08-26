#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Frankie
@date: 20180308
"""
from ..core import Step, Configure
from .Cellranger import Cellranger
import os

class Seuratrun(Step):
    """
        Seuratrun is a Step to tackle with Seurat object and finish Seurat analysis.
        See __init__ to initialize this step.
        > Seuratrun(): __init__  parameters
            inputdir: str or str list
                a directory contain seurat object, for example: /home/data/Seuratpreprocessing
            outputdir: str
                outputdir path of all the results, default: ~/step_02_Seuratrun/
            rscript: str
                rscript file path, for example: /home/data/Seuratrun.R
            xlowcut: num
                low threshold of variable genes mean expression
            xhighcut: num
                high threshold of variable genes mean expression
            ycut: num
                threshold of variable genes dispersion
            cmdParam: str or list of string
                current unsupported
         > Seuratrun()(): __call__ parameters
            Available upstream objects combinations:
            (Seuratpreprocessing)
            (to be added)
        """
    def __init__(self,
                 outputdir=None,
                 # barcodes=None,
                 # genes=None,
                 inputdir = None,
                 rscript=None,
                 xlowcut = None,
                 xhighcut = None,
                 ycut = None,
                 # datatype = None,run
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
        # self.setParam('datatype', datatype)
        self.setParam('xlowcut', xlowcut)
        self.setParam('xhighcut', xhighcut)
        self.setParam('ycut', ycut)

        self.initIO()
        self._setMultiRun()

    def impInitIO(self, ):
        # matrix = self.getParamIO('matrix')
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        outputdir = self.getParamIO('outputdir')
        inputdir = self.getParamIO('inputdir')

        self.setInputRscript('rscript', '../rscript/Seuratrun.R')
        # if outputdir is None, os will error
        if outputdir is None:
           self.setParamIO('outputdir', Configure.getTmpDir())
           outputdir = self.getParamIO('outputdir')

        # set output/input paths
        # if inputdir is None, os will error
        if inputdir is not None:
            self.setInputDirOrFile('inputfile', os.path.join(inputdir, 'Preprocessing.Rdata'))
            # self.setInputDirOrFile('barcodes', os.path.join(inputdir, 'barcodes.tsv'))
            # self.setInputDirOrFile('genes', os.path.join(inputdir, 'genes.tsv'))
            # self.setInputDirOrFile('matrix', os.path.join(inputdir, 'matrix.mtx'))
            # self.setOutputDirNTo1('violinplot', os.path.join(outputdir, 'violinplot.jpeg'), '', 'barcodes')
            # self.setOutputDirNTo1('geneplot', os.path.join(outputdir, 'geneplot.jpeg'), '', 'barcodes')
            self.setOutputDirNTo1('Elbowplot', os.path.join(outputdir, 'Elbowplot.jpeg'), '', 'inputfile')
            self.setOutputDirNTo1('TSNEplot', os.path.join(outputdir, 'TSNEplot.jpeg'), '', 'inputfile')

    def call(self, *args):
        # the first object
        cellrangerUpstream = args[0]
        # set all required input parameters from upstream objects
        # if isinstance(cellrangerUpstream, Cellranger):
        #     datatype = '10x'
        #     self.setParam('datatype', datatype)
        self.setParamIO('inputdir', cellrangerUpstream.getParamIO('outputdir'))
        # self.setParamIO('genes', cellrangerUpstream.getOutput('genes'))
        # self.setParamIO('matrix', cellrangerUpstream.getOutput('matrix'))

        # print(self.cmdline)

    def _multiRun(self, ):
        # get input parameters
        # barcodes = self.getParamIO('barcodes')
        # genes = self.getParamIO('genes')
        # matrix = self.getParamIO('matrix')
        inputfile = self.getInput('inputfile')
        rscript = self.getInput('rscript')
        # datatype = self.getParam('datatype')
        # get output parameters
        outputdir = self.getParamIO('outputdir')

        xlow = self.getParam('xlowcut')
        xhigh = self.getParam('xhighcut')
        y = self.getParam('ycut')

        if (xlow is None) and (xhigh is None) and (y is None):
            cmdline = ['Rscript',
                       rscript,
                       inputfile,
                       # datatype,
                       outputdir, 'xlow', 'xhigh', 'y']
            self.callCmdline(cmdline=cmdline, dockerVersion='V1')
        else:
            cmdline = ['Rscript',
                       rscript,
                       inputfile,
                       # datatype,
                       outputdir, xlow, xhigh, y]
            self.callCmdline(cmdline=cmdline, dockerVersion='V1')
        #
        # replace file path
        # inputdir = inputdir.replace('/home/cfeng', '/data')
        # # matrix = matrix.replace('/home/cfeng', '/data')
        # rscript = rscript.replace('/home/cfeng', '/data')
        # # barcodes = barcodes.replace('/home/cfeng', '/data')
        #
        # outputdir = outputdir.replace('/home/cfeng', '/data')
        # if datatype is ('10x' or '10X' or 'tenx' or 'Tenx' or 'TenX' or 'tenX'):
        #     cmdline = ['Rscript',
        #                rscript,
        #                inputdir,
        #                #datatype,
        #                outputdir]
        #     print(cmdline)
        #     self.callCmdline(cmdline=cmdline, dockerVersion='V1')
        # else:
        #     cmdline = ['Rscript',
        #                rscript,
        #                inputdir,
        #                # datatype,
        #                outputdir]
        #     print(cmdline)
        #     self.callCmdline(cmdline=cmdline, dockerVersion='V1')

    def getMarkdownEN(self,):
        rmd = '''
## Seurat analysis
`Runobject = Seuratrun(inputdir, outputdir, rscript)(upstream)`


### Determine statistically significant principal components

To overcome the extensive technical noise in any single gene for scRNA-seq data, determining how many PCs to include downstream is therefore an important step.

![]({elbowimagepath})

### Run Non-linear dimensional reduction (tSNE)

Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets.

![]({tsneimagepath})
        '''.format(elbowimagepath=self.getOutput('Elbowplot'),
                   tsneimagepath=self.getOutput('TSNEplot'))
        return rmd
