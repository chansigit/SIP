# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/10 8:49
@Author  : Weizhang
@FileName: SRAToFastq1.py
"""

from ..core import Step, Configure


class SRAToFastq(Step):
    def __init__(self,
                 sraInput=None,
                 fastqOutputDir=None,
                 cmdParam=None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, **kwargs)

        # set IO parameters
        self.setParamIO('sraInput', sraInput)
        self.setParamIO('fastqOutputDir', fastqOutputDir)

        self.initIO()

        # set other parameters
        self.setParam('fastqc', True)

    def impInitIO(self):
        sraInput = self.getParamIO('sraInput')
        fastqOutputDir = self.getParamIO('fastqOutputDir')

        # set all input files
        self.setInputDirOrFile('sraInput', sraInput)
        # set all output files
        self.setOutputDir1To1('fastqOutput1', fastqOutputDir, None, '_1.fastq', 'sraInput', '')
        self.setOutputDir1To1('fastqOutput1_fqch', fastqOutputDir, None, '_1_fastqc.html', 'sraInput', '')
        self.setOutputDir1To1('fastqOutput1_fqczip', fastqOutputDir, None, '_1_fastqc.zip', 'sraInput', '')
        self.setOutputDir1To1('fastqOutput2', fastqOutputDir, None, '_2.fastq', 'sraInput', '')
        self.setOutputDir1To1('fastqOutput2_fqch', fastqOutputDir, None, '_2_fastqc.html', 'sraInput', '')
        self.setOutputDir1To1('fastqOutput2_fqczip', fastqOutputDir, None, '_2_fastqc.zip', 'sraInput', '')

        if fastqOutputDir is None:
            self.setParamIO('fastqOutputDir', Configure.getTmpDir())

        if sraInput is not None:
            self._setInputSize(len(self.getInputList('sraInput')))

    def call(self, *args):

        print("SRAToFastq has no upstream!!!\n")

    def _multiRun(self,):
        pass

    def _singleRun(self, i):
        sraInput = self.getInputList('sraInput')
        fastqOutputDir = self.getParamIO('fastqOutputDir')
        fastqOutput1 = self.getOutputList('fastqOutput1')
        fastqOutput1_fqch = self.getOutputList('fastqOutput1_fqch')
        fastqOutput1_fqczip = self.getOutputList('fastqOutput1_fqczip')
        fastqOutput2 = self.getOutputList('fastqOutput2')
        fastqOutput2_fqch = self.getOutputList('fastqOutput2_fqch')
        fastqOutput2_fqczip = self.getOutputList('fastqOutput2_fqczip')

        if self.getParam('fastqc'):  # do fastqc
            cmdline1 = ['fastq-dump', '--split-3',
                        sraInput[i],
                        '-O', fastqOutputDir
                        ]
            result = self.callCmdline('V1', cmdline1)

            cmdline2 = ['fastqc', fastqOutput1[i], fastqOutput2[i]]
            result = self.callCmdline('V1', cmdline2)
        else:  # do not do fastqc
            cmdline1 = ['fastq-dump', '--split-3',
                        sraInput[i],
                        '-O', fastqOutputDir
                        ]
            result = self.callCmdline('V1', cmdline1)


    def getMarkdownEN(self,):
        fastqOutput1_fqch = self.getOutputList('fastqOutput1_fqch')
        fastqOutput2_fqch = self.getOutputList('fastqOutput2_fqch')

        fastqc1_path = []
        for qc in fastqOutput1_fqch:
            fastqc1_path.append("\"" + qc + "\"")
        fastqc1_path_str = ",".join([i for i in fastqc1_path])
        fastqc1_path_str = "c(" + fastqc1_path_str + ")"

        fastqc2_path = []
        for qc in fastqOutput2_fqch:
            fastqc2_path.append("\"" + qc + "\"")
        fastqc2_path_str = ",".join([i for i in fastqc2_path])
        fastqc2_path_str = "c(" + fastqc2_path_str + ")"

        mdtext ="""
## SRAToFastq and FASTQC Result
```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
library(knitr)
library(kableExtra)
fq1 <- {fq1path}
fq2 <- {fq2path}
fq <- cbind(fq1, fq2)
colnames(fq) <- c("Quality Control for _1.fastq", "Quality Control for _2.fastq")
```
The FASTQC html report path are as follows:

```{{r eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE}}
kable(fq, "html") %>% kable_styling() %>% scroll_box(width = "1100px", height = "500px")
```


        """.format(fq1path=fastqc1_path_str, fq2path=fastqc2_path_str)
        return mdtext
