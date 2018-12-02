from ..core import Step, Configure
import subprocess
import os

class StarAlign(Step):
    def __init__(self,
                 fastqInput = None,
                 outFileDir = None,
                 genomeDir = None,
                 outFileNamePrefix = None,
                 threads = None,
                 #outSamType = 'BAM'
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam, ** kwargs)

        Configure.setRefSuffix('starRef', '_refdata-cellranger-2.1.0/star', check=False)
        self.setParamIO('fastqInput', fastqInput)
        self.setParamIO('outFileDir', outFileDir)
        self.setParamIO('genomeDir', genomeDir)

        bamOutput = None
        logFinalOutput = None
        logProgressOutput = None
        logOutput = None
        tabOutput = None
        if outFileDir is not None:
            bamOutput = os.path.join(outFileDir, outFileNamePrefix+'Aligned.out.sam')
            logFinalOutput = os.path.join(outFileDir, outFileNamePrefix+'Log.final.out')
            logProgressOutput = os.path.join(outFileDir, outFileNamePrefix+'Log.progress.out')
            logOutput = os.path.join(outFileDir, outFileNamePrefix+'Log.out')
            tabOutput = os.path.join(outFileDir, outFileNamePrefix+'tabOutput')
        self.setParamIO('bamOutput', bamOutput)
        self.setParamIO('logFinalOutput', logFinalOutput)
        self.setParamIO('logProgressOutput', logProgressOutput )
        self.setParamIO('logOutput', logOutput)
        self.setParamIO('tabOutput', tabOutput)

        if threads is None:
            threads = Configure.getThreads()
        self.setParam('threads', threads)
        self.setParam('outFileNamePrefix', outFileNamePrefix)
        #self.setParam('outSamType', outSamType)

        self.initIO()

    def impInitIO(self,):
        fastqInput = self.getParamIO('fastqInput')
        outFileDir = self.getParamIO('outFileDir')
        bamOutput = self.getParamIO('bamOutput')
        logFinalOutput = self.getParamIO('logFinalOutput')
        logProgressOutput = self.getParamIO('logProgressOutput')
        logOutput = self.getParamIO('logOutput')
        tabOutput = self.getParamIO('tabOutput')
        genomeDir = self.getParamIO('genomeDir')
        outFileNamePrefix = self.getParam('outFileNamePrefix')

        if genomeDir is None:
            self.setParamIO('genomeDir', Configure.getConfig('starRef'))
            genomeDir = self.getParamIO('genomeDir')

        for i in ['chrLength.txt', 'chrName.txt',
                  'exonGeTrInfo.tab', 'geneInfo.tab',
                  'genomeParameters.txt', 'SAindex',
                  'sjdbList.fromGTF.out.tab', 'transcriptInfo.tab',
                  'chrNameLength.txt', 'chrStart.txt',
                  'exonInfo.tab', 'Genome', 'SA',
                  'sjdbInfo.txt', 'sjdbList.out.tab']:
            self.setInputDirOrFile(i, os.path.join(genomeDir, i))

        self.setInputDirOrFile('fastqInput', fastqInput)
        self.setOutputDirNTo1('bamOutput', bamOutput, outFileNamePrefix+'Aligned.out.sam', 'fastqInput')
        self.setOutputDirNTo1('logFinalOutput', logFinalOutput, outFileNamePrefix+'Log.final.out', 'fastqInput')
        self.setOutputDirNTo1('logProgressOutput', logProgressOutput, outFileNamePrefix+'Log.progress.out', 'fastqInput')
        self.setOutputDirNTo1('logOutput', logOutput, outFileNamePrefix+'Log.out', 'fastqInput')
        self.setOutputDirNTo1('tabOutput', tabOutput, outFileNamePrefix+'SJ.out.tab', 'fastqInput')

        if fastqInput is not None:
            self._setInputSize(len(self.getInputList('fastqInput')))

        if outFileDir is None:
            self.setParamIO('outFileDir', Configure.getTmpDir())
        if bamOutput is None:
            self.setParamIO('bamOutput', Configure.getTmpPath(outFileNamePrefix+'Aligned.out.sam'))
        if logFinalOutput is None:
            self.setParamIO('logFinalOutput', Configure.getTmpPath(outFileNamePrefix+'Log.final.out'))
        if logProgressOutput is None:
            self.setParamIO('logProgressOutput', Configure.getTmpPath(outFileNamePrefix+'Log.progress.out'))
        if logOutput is None:
            self.setParamIO('logOutput', Configure.getTmpPath(outFileNamePrefix+'Log.out'))
        if tabOutput is None:
            self.setParamIO('tabOutput', Configure.getTmpPath(outFileNamePrefix+'SJ.out.tab'))

    def call(self, *args):
        fastqUpstream = args[0]
        self.setParamIO('fastqInput', fastqUpstream.getOutput('fastqOutput'))

    def _singleRun(self, i):
        fastqInput = self.getInputList('fastqInput')
        outFileDir = self.getParamIO('outFileDir')
        genomeDir = self.getParamIO('genomeDir')

        threads = self.getParam('threads')
        outFileNamePrefix = self.getParam('outFileNamePrefix')
        cmdline = [
                 '/root/software/STAR/source/STAR', '--runThreadN %d'%(threads),'--genomeDir %s'%(genomeDir),
                '--readFilesIn %s'%(fastqInput[i]), '--outFileNamePrefix %s'%(os.path.join(outFileDir, outFileNamePrefix))
        ]
        self.callCmdline('V1', cmdline)

    def getMarkdownEN(self,):
        mdtext="""
## StarAlign Result
The alignment result of STAR is shown below:
```{{r echo=FALSE}}
slfo <- file("{logFinalOutput}", "r", blocking=FALSE)
lines <- readLines(slfo)
strlist <- strsplit(lines, split='\t')
strmtx <- do.call(rbind, strlist)
strmtx[c(7,22,27,31), 2] <- ""
strmtx <- rbind(strmtx[1:4,], c("", ""), strmtx[5:33,])
data.frame(Info=strmtx[,1], Data=strmtx[,2])
```

        """.format(logFinalOutput=self.getOutput('logFinalOutput'))
        return mdtext
