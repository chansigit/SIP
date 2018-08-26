# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 19:45:01 2018

@author: WeiZheng
"""

import os
import subprocess
from .StepBase import Configure, Schedule, Report

class Flow:
    def __init__(self, resultDir='.', refdir=None, genome=None, threads=None):
        self.refdir = refdir
        self.genome = genome
        self.threads = threads
        if resultDir is None:
            raise Exception('resultDir can not be None')
        self.resultDir = os.path.abspath(resultDir)
        
        self.objs = dict()
        self.param = dict()
        
    
    def __call__(self,*args):
        #if len(args) == 0:
        #    raise Exception('there should be one upsteram Pipe at least')
        for i in range(len(args)):
            if not isinstance(args[i],Flow):
                raise Exception('only Pipe subclasses are supported')
        if len(args) > 0 :
            self._call(*args)
        self.run()
        return self
    
    def _call(self, *args):
        raise Exception("method: call must be overwrited")
    
    def run(self,):
        if self.refdir is not None:
            Configure.setRefDir(self.refdir)
        if self.genome is not None:
            Configure.setGenome(self.genome)
        if self.threads is not None:
            Configure.setThreads(self.threads)
            
        self.tmpdir = os.path.join(self.resultDir,'intermediate_results')
        self.reportdir = os.path.join(self.resultDir,'report')
        self.finalResultDir = os.path.join(self.resultDir,'final_results')
        os.makedirs(self.resultDir, exist_ok=True)
        os.makedirs(self.tmpdir, exist_ok=True)
        os.makedirs(self.reportdir, exist_ok=True)
        os.makedirs(self.finalResultDir, exist_ok=True)
        tmpdir_store = Configure.getTmpDir()
        Configure.setTmpDir(self.tmpdir)
        self._build()
        Schedule.run(report = False)
        #Schedule.stopDocker()
        Configure.setTmpDir(tmpdir_store) 
        count = 0
        for key in self.objs:
            obj = self.objs[key]
            if isinstance(obj,Report):
                print(obj)
                self._linkRecursive(obj.getOutput('reportHTML'),
                                    os.path.join(self.getReportDir(),'report'+str(count)+'.html'))
                count +=1
                self._linkRecursive(os.path.join(os.path.dirname(obj.getOutput('reportHTML')),'links'),
                                    os.path.join(self.getReportDir(),'links'))
                print(count)
        self._copy()
        
        
    def _build(self,):
        raise Exception('_build(self,) must be implemented')
    
    def _copy(self,):
        raise Exception('_copy(self,) must be implemented')

    def _setObj(self, objName, obj):
        self.objs[objName] = obj
    
    def _getObj(self, objName):
        return self.objs[objName]
    
    def _setParam(self, paramName,value):
        self.param[paramName] = value
    
    def _getParam(self, paramName):
        return self.param[paramName]
    
    def getReportDir(self,):
        return self.reportdir
    
    def getMedRsDir(self,):
        return self.tmpdir
    
    def getFinalRsDir(self,):
        return self.finalResultDir
    
    def getFinalRsPath(self, fileOrFolderName):
        return os.path.join(self.finalResultDir,fileOrFolderName)
    
    def _linkRecursive(self,srcPath,desPath,desFolder=True):
        #print(originPath)
        #print('To')
        #print(desPath)
        if not desPath.startswith(self.finalResultDir):
            desPath = self.getFinalRsPath(desPath)
        originPaths = srcPath
        if isinstance(srcPath,list):
            os.makedirs(desPath, exist_ok=True)
        else:
            if desFolder:
                os.makedirs(desPath, exist_ok=True)
            originPaths = [srcPath]
        for originPath in originPaths:
            print('------------')
            print(originPath)
            print('------------')
            if os.path.isfile(originPath):
                #subprocess.run(['rm','-rf',desPath])
                subprocess.run(['ln','-f',originPath,desPath],check=True)
                continue
            elif os.path.isdir(originPath):
                os.makedirs(desPath,exist_ok=True)
                files = os.listdir(originPath)
                for f in files:
                    self._linkRecursive(os.path.join(originPath,f),os.path.join(desPath,f),desFolder=False)
            else:
                #print(['originPath:',originPath,'is neither dir nor file'])
                raise Exception('originPath:',originPath,'is neither dir nor file')
