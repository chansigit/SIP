# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 11:12:29 2018

@author: WeiZheng
"""
from ..core import Flow, Report, Configure, Schedule
from ..steps import AdapterRemoval, Bowtie2 

class FlowExample(Flow):
    def __init__(self,
                 fastqInput1,
                 fastqInput2,
                 refdir, 
                 genome, 
                 resultDir='./result',                  
                 threads=None,                 
                 adapter1 = None,
                 adapter2 = None,
                 X = 2000,):
        super(FlowExample,self).__init__(resultDir=resultDir, 
                                         refdir=refdir, 
                                         genome=genome, 
                                         threads=threads)
        self._setParam('fastqInput1',fastqInput1)
        self._setParam('fastqInput2',fastqInput2)
        self._setParam('adapter1',adapter1)
        self._setParam('adapter2',adapter2)
        self._setParam('X',X)
        
    def _call(self,*args):
        # args[0]._getObj('FastqDump').getOutput('fastqOutput')
        pass
    
    def _build(self,):
        fastqInput1 = self._getParam('fastqInput1')
        fastqInput2 = self._getParam('fastqInput2')
        adapter1 = self._getParam('adapter1')
        adapter2 = self._getParam('adapter2')
        X = self._getParam('X')
        
        adrm = AdapterRemoval(fastqInput1=fastqInput1, fastqInput2=fastqInput2,
                              adapter1=adapter1, adapter2=adapter2)
        bt = Bowtie2(X=X)(adrm)
        
        rp = Report()
        rp.add('Section for AdapterRemoval',[adrm])
        rp.add('Section for Bowtie2',[bt])
        
        self._setObj('AdapterRemoval',adrm)
        self._setObj('Mapping',bt)
        self._setObj('Report',rp)
        
    def _copy(self,):
        
        self._linkRecursive(self._getObj('Mapping').getOutput('samOutput'),
                            self.getFinalRsPath('MappedReads'))
