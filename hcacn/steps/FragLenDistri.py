# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/30 17:53
@Author  : Weizhang
@FileName: FragLenDistri.py

generate fragment length distribution plot for scATAC-seq

"""


from ..core import Step, Configure
import os.path

class FragLenDistri(Step):
    def __init__(self,
                 bedInput=None,
                 figureOutputDir=None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('bedInput', bedInput)
        self.setParamIO('figureOutputDir', figureOutputDir)
        self.initIO()

        # set other parameters

    def impInitIO(self):
        bedInput = self.getParamIO('bedInput')
        figureOutputDir = self.getParamIO('figureOutputDir')

        # FragLenDistri.R
        self.setInputRscript('FragLenDistriR', 'FragLenDistri.R')

        # set all input files
        self.setInputDirOrFile('bedInput', bedInput)
        # set all output files
        self.setOutputDir1To1('figureOutput', figureOutputDir, None, 'bmp', 'bedInput')

        if figureOutputDir is None:
            self.setParamIO('figureOutputDir', Configure.getTmpDir())

        if bedInput is not None:
            self._setInputSize(len(self.getInputList('bedInput')))

    def call(self, *args):
        bedUpstream = args[0]
        self.setParamIO('bedInput', bedUpstream.getOutput('bedOutput'))

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        bedInput = self.getInputList('bedInput')  # list
        figureOutput = self.getOutputList('figureOutput')  # a file name
        rScript = self.getInput('FragLenDistriR')

        cmdline = [
            'Rscript', rScript,
            bedInput[i], figureOutput[i]
        ]
        result = self.callCmdline('V1', cmdline)

    def getMarkdownEN(self, ):

        figureOutput = self.getOutputList('figureOutput')

        figure_path = []
        for fi in figureOutput:
            figure_path.append("\"" + fi + "\"")
        figure_path_str = ",".join([i for i in figure_path])
        fastqc1_path_str = "c(" + figure_path_str + ")"

        mdtext = """
## Fragment Length Distribution Result

The following is an example about fragment length distribution for scATAC-seq data.

```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
library(knitr)
library(kableExtra)
figure_path <- {figureOutput}
example_figure <- figure_path[1]

figure_path <- as.data.frame(figure_path)
colnames(figure_path) <- "figure file path"
```

![](`r example_figure`)

For more information, please seek in the following table:

```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
kable(figure_path, "html") %>% kable_styling() %>% scroll_box(width = "800px", height = "500px")
```

        """.format(figureOutput=fastqc1_path_str)
        return mdtext
