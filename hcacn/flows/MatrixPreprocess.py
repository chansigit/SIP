# -*- coding: utf-8 -*-
"""
@author: Zhenyi Wang
2018/3/29
"""
from ..core import Flow,Report,Configure,Schedule
from ..steps import matrixRw

class MatrixPreprocess(Flow):
    def __init__(self,
                 matrixdata,
                 outputpath,
                 resultDir='./result'):
        super(MatrixPreprocess,self).__init__(resultDir=resultDir, 
                                         refdir=None, 
                                         genome=None, 
                                         threads=None)

        self._setParam('matrixdata',matrixdata)
        self._setParam('outputpath',outputpath)
       

    def _call(self,*args):
        # args[0]._getObj('FastqDump').getOutput('fastqOutput')
        pass
    
    def _build(self,):
        matrixdata = self._getParam('matrixdata')
        outputpath = self._getParam('outputpath')
    
        mRW = matrixRw(matrixdata = matrixdata, outputpath = outputpath)       
      
        rp = Report()
        rp.add('Section for Matrix pre-process',[mRW])
        self._setObj('MatrixPreprocess',mRW)
        self._setObj('Report',rp)
        
    def _copy(self,):
        self._linkRecursive(self._getObj('MatrixPreprocess').getOutput('processedMatrix'),
                            self.getFinalRsPath('processedMatrix'))
        
