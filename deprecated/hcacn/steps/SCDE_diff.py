from ..core import Step,Configure
import os

class SCDE_fidd(Step):
	def __init__(self,
		matrixdata = None,
		outputpath = None,
		cmdParam = None,
		**kwargs):

	super(Step, self).__init__(cmdParam,**kwargs)
        # set all input and output parameters
        self.setParamIO('matrixdata',sceInput)
        self.setParamIO('outputpath',outputpath)
        #self.setParamIO('clusterName1',clusterName1)
        #self.setParamIO('clusterName2',clusterName2)
        # call self.initIO()
        self.initIO()
        #set other parameters
	def impInitIO(self,):
        
        matrixdata = self.getParamIO('matrixdata')
        outputpath = self.getParamIO('outputpath')
        #clusterName1 = self.getParamIO('clusterName1')
        #clusterName2 = self.getParamIO('clusterName2')
        if outputpath is None:
        	self.setParamIO('outputpath',Configure.getTmpDir())
        	outputpath = self.getParamIO('outputpath')  

        self.setInput('matrixdata',matrixdata)

        self.setOutputDir1To1('scdeDiff',outputpath,'scdeDiff','txt','matrixdata')

    def call(self, *args):
    	pass

    	def _multiRun(self,):
        matrixdata = self.getInput('matrixdata')
        cmdline = ['Rscript SCDE_diff.R',
        			matrixdata,
        			self.getParamIO('outputpath')
                    #clusterName1,
                    #clusterName2
        ]
        print(''.join(cmdline))
        self.callCmdline('V1',cmdline)