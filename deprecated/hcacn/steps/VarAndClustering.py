# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 17:13
@Author  : Weizhang
@FileName: VarAndClustering.py
"""


from ..core import Step, Configure
from .GenPeakWithFilter import GenPeakWithFilter
import os.path


class VarAndClustering(Step):
    def __init__(self,
                 bamInput=None,
                 figureOutput=None,
                 peakInput=None,
                 genome=None,
                 threads=1,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bamInput', bamInput)
        self.setParamIO('figureOutput', figureOutput)
        self.setParamIO('peakInput', peakInput)

        self.initIO()

        # set other parameters
        self.setParam('threads', threads)
        if genome is None:
            self.setParam('genome', Configure.getGenome())
        else:
            self.setParam('genome', genome)

        self._setUpstreamSize(2)

    def impInitIO(self):
        bamInput = self.getParamIO('bamInput')
        figureOutput = self.getParamIO('figureOutput')
        peakInput = self.getParamIO('peakInput')

        # ChromVarUtility.R
        self.setInputRscript('ChromVarUtilityR', 'ChromVarUtility.R')


        # set all input files
        self.setInputDirOrFile('bamInput', bamInput)
        self.setInputDirOrFile('peakInput', peakInput)
        # set all output files
        self.setOutputDirNTo1('var.tiff', None, 'variation.bmp', 'bamInput')
        self.setOutputDirNTo1('clustring.tiff', None, 'clustring.bmp', 'bamInput')
        self.setOutputDirNTo1('dev.matrix', None, 'devmatrix.txt', 'bamInput')
        self.setOutputDirNTo1('var.matrix', None, 'varmatrix.txt', 'bamInput')
        self.setOutputDirNTo1('clu.matrix', None, 'clumatrix.txt', 'bamInput')

        if peakInput is not None:
            self._setInputSize(len(self.getInputList('peakInput')))

    def call(self, *args):
        bamUpstream = args[0]
        peakUpstream = args[1]

        self.setParamIO('bamInput', bamUpstream.getOutputList('bamOutput'))
        self.setParamIO('peakInput', peakUpstream.getOutputList('bedOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        bamInput = self.getInputList('bamInput')
        peakInput = self.getInputList('peakInput')
        rScript = self.getInput('ChromVarUtilityR')
        varOutput = self.getOutputList('var.tiff')
        clustringOutput = self.getOutputList('clustring.tiff')
        devpath = self.getOutputList('dev.matrix')
        varpath = self.getOutputList('var.matrix')
        clupath = self.getOutputList('clu.matrix')
        cmdline = [
            'Rscript', rScript,
            os.path.dirname(bamInput[i]),
            peakInput[i],
            str(self.getParam('threads')),
            str(self.getParam('genome')),
            varOutput[i], clustringOutput[i], devpath[i], varpath[i], clupath[i]
        ]
        result = self.callCmdline('V1', cmdline)

    def getMarkdownEN(self, ):
        varOutput = self.getOutputList('var.tiff')
        clustringOutput = self.getOutputList('clustring.tiff')
        devpath = self.getOutputList('dev.matrix')
        varpath = self.getOutputList('var.matrix')
        clupath = self.getOutputList('clu.matrix')

        mdtext = """
## Variation and Clustering Result

The following is variation and clustering plots.
```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
var_path <- "{varOutput}"
clu_path <- "{clustringOutput}"
devpath <- "{devpath}"
varpath <- "{varpath}"
clupath <- "{clupath}"
```

![Variation Plot](`r var_path`)

![Clustering Plot](`r clu_path`)

You can inspect the data from the following path:

For deviation matrix: `r devpath`

For variability matrix: `r varpath`

For clustring matrix: `r clupath`

                """.format(varOutput=varOutput[0], clustringOutput=clustringOutput[0],
                           devpath=devpath[0], varpath=varpath[0], clupath=clupath[0])
        return mdtext
