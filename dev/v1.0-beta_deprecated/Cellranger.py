
from ..core import Step, Configure
import os

class Cellranger(Step):
    """
        Cellranger is a Step to count 10X raw data to UMI counts.
        See __init__ to initialize this step.
        > Cellranger(): __init__  parameters
            fastqInput: str or str list
                a directory contain all of fastq files, for example: /home/data/fastqs
            outputdir: str
                outputdir path of all the cellranger results, default: ~/step_00_Cellranger/
            refile: str
                cellranger transcriptome file path, for example: /home/data/refdata-cellranger-hg19_and_mm10-1.2.0
            expectcells: int
                the expect number of cells
            cmdParam: str or list of string
                current unsupported
        """

    def __init__(self,
                 outputdir = None,
                 fastqInput = None,
                 refile = None,
                 expectcells = None,
                 cmdParam = None,
                 **kwargs):
        super(Step, self).__init__(cmdParam,**kwargs)
        Configure.setRefSuffix('cellranger', '.refdata-cellranger-hg19-1.2.0', check=False)
        self.setParamIO('fastqInput', fastqInput)
        self.setParamIO('refile', refile)
        self.setParamIO('outputdir', outputdir)
        self.initIO()

        self.setParam('expectcells', expectcells)
        #Configure.enableDocker(False)
        self._setMultiRun()

    def impInitIO(self,):
        fastqInput = self.getParamIO('fastqInput')
        refile = self.getParamIO('refile')
        outputdir = self.getParamIO('outputdir')

        self.setInputDirOrFile('fastqInput', fastqInput)

        if refile is None:
            self.setParamIO('refile', Configure.getConfig('cellranger'))
            refile = self.getParamIO('refile')

        self.setInputDirOrFile('version', os.path.join(refile, 'version'))
        self.setInputDirOrFile('Reference', os.path.join(refile, 'reference.json'))
        self.setInputDirOrFile('README', os.path.join(refile, 'README.BEFORE.MODIFYING'))
        for i in ['chrLength.txt', 'chrName.txt',
                  'exonGeTrInfo.tab', 'geneInfo.tab',
                  'genomeParameters.txt', 'SAindex',
                  'sjdbList.fromGTF.out.tab', 'transcriptInfo.tab',
                  'chrNameLength.txt', 'chrStart.txt',
                  'exonInfo.tab', 'Genome', 'SA',
                  'sjdbInfo.txt', 'sjdbList.out.tab']:
            self.setInputDirOrFile(i, os.path.join(refile, 'star', i))
        self.setInputDirOrFile('genes.pickle', os.path.join(refile, 'pickle', 'genes.pickle'))
        self.setInputDirOrFile('genes.gtf', os.path.join(refile, 'genes', 'genes.gtf'))
        self.setInputDirOrFile('genome.fa', os.path.join(refile, 'fasta', 'genome.fa'))

        if outputdir is None:
            self.setParamIO('outputdir', Configure.getTmpDir())
            outputdir = self.getParamIO('outputdir')
            self.resultdir = 'Cellranger'
        else:
            self.resultdir = ''

        self.setParamIO('finaldir', os.path.join(outputdir, self.resultdir, 'outs', 'filtered_gene_bc_matrices', 'hg19'))
        self.setOutputDirNTo1('summary', os.path.join(outputdir, self.resultdir,'outs', 'web_summary.html'), '', 'fastqInput')
        self.setOutputDirNTo1('genes', os.path.join(outputdir, self.resultdir,'outs', 'filtered_gene_bc_matrices', 'hg19', 'genes.tsv'), '', 'fastqInput')
        self.setOutputDirNTo1('matrix', os.path.join(outputdir, self.resultdir,'outs', 'filtered_gene_bc_matrices', 'hg19', 'matrix.mtx'), '', 'fastqInput')
        self.setOutputDirNTo1('barcodes', os.path.join(outputdir,self.resultdir, 'outs', 'filtered_gene_bc_matrices', 'hg19', 'barcodes.tsv'), '', 'fastqInput')

    def call(self, *args):
        pass
        # print(self.cmdline)

    def _multiRun(self, ):
        fastqInput = self.getParamIO('fastqInput')
        refile = self.getParamIO('refile')
        # id = self.getParam('id')
        outputdir = self.getParamIO('outputdir')
        expectcells = self.getParam('expectcells')
        # cd outputdir to run the command
        # self.cmdline1 = ['cd', '%s' % outputdir]
        # self.callCmdline(self.cmdline1)
        if expectcells is not None:
            if self.resultdir is '':
                # use given path
                dir = outputdir.split('/')[-1]
                if os.path.isdir(outputdir):
                    cmdline = ['rm -r', '%s' % outputdir]
                    self.callCmdline(cmdline=cmdline, dockerVersion='V1',stdoutToLog=False)
                cmdline = ['cellranger', 'count','--id=%s --expect-cells=%s --transcriptome=%s --fastqs=%s'
                           % (dir, str(expectcells), refile, fastqInput)]
                print(cmdline)
                self.callCmdline(cmdline= cmdline, dockerVersion='V1', stdoutToLog=False)
            else:
                # use default filepath
                cmdline = ['cd', outputdir, '&&', 'rm -rf Cellranger', '&&', 'cellranger', 'count',
                           '--id=Cellranger', '--expect-cells=%s --transcriptome=%s --fastqs=%s'
                           % (str(expectcells), refile, fastqInput)]
                self.callCmdline(cmdline=cmdline, dockerVersion='V1', stdoutToLog=False)

        else:
            if self.resultdir is '':
                # use given path
                dir = outputdir.split('/')[-1]
                if os.path.isdir(outputdir):
                    cmdline = ['rm -r', '%s' % outputdir]
                    self.callCmdline(cmdline=cmdline, stdoutToLog=False)
                cmdline = ['cellranger', 'count', '--id=%s  --transcriptome=%s --fastqs=%s'
                           % (dir, refile, fastqInput)]
                print(cmdline)
                self.callCmdline(cmdline=cmdline, dockerVersion='V1', stdoutToLog=False)
            else:
                # use default filepath
                # cd target dir
                # run 10x cellranger
                cmdline = ['cd', outputdir,'&&','rm -rf Cellranger', '&&', 'cellranger', 'count',
                           '--id=Cellranger', '--transcriptome=%s --fastqs=%s'
                           % (refile, fastqInput)]
                self.callCmdline(cmdline=cmdline, dockerVersion='V1', stdoutToLog=False)


    def getMarkdownEN(self, ):
        rmd = '''	
## Cellranger

`Cellrangerobject = Cellranger(fastqinput, outputdir, refile, expectcells)`

### Cellranger count fastqs to UMIs

- Final outputs will be saved in ~/step_XX_Cellranger/outs/
- Matrix file in sparse format can be found in /outs/filtered_gene_bc_matrices/hg19
- Summary.html contain most of data infos can be found in /outs/web_summary.html 


        '''
        return rmd


