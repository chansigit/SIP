# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 19:22:29 2018
@author: WeiZheng
"""

import numpy as np

class ValidGraph:
    def __init__(self,edges,nodes):
        self.nodeIdx = self.genNodeIdx(edges, nodes)
        self.nodeAttr = {}        
        self.adjMatrix = self.getAdjMatrix(self.nodeIdx,edges)
        
    def genNodeIdx(self,edges, nodes):
        nodeIdx = {}
        idx = 0
        for i in range(len(edges)):
            if len(edges[i]) != 2:
                raise Exception("edges",i,"do not contain 2 node")
            for j in range(2):
                if not isinstance(edges[i][j],str):
                    raise Exception("node",edges[i][j],"should be an string")
                if not edges[i][j] in nodeIdx:
                    nodeIdx[edges[i][j]] = idx
                    idx += 1
        for node in nodes:
            if not isinstance(node,str):
                raise Exception("node",node,"should be an string")
            if node not in nodeIdx:
                nodeIdx[node] = idx
                idx +=1
                    
        return nodeIdx
    
    def getAdjMatrix(self,nodeIdx,edges):
        nodeNum = len(nodeIdx)
        adjMatrix = np.zeros((nodeNum,nodeNum),dtype=bool)
        for i in range(len(edges)):
            adjMatrix[nodeIdx[edges[i][0]],nodeIdx[edges[i][1]]] = True
        return adjMatrix
        
    def setNodeAttrs(self,attrType,aDict):
        attrDict = self.nodeIdx.copy()
        for key in aDict.keys():
            attrDict[key] = aDict[key]
        self.nodeAttr[attrType] =  attrDict
    
    def isConnect(self,fromNode, toNode): 
        return self.adjMatrix[self.getNodeIdx(fromNode),self.getNodeIdx(toNode)]
    
    def getNodeIdx(self,nodeStr):
        return self.nodeIdx[nodeStr] 
    
    def getNodeNum(self,):
        return self.adjMatrix.shape[0]
     


class GraphMng:
    def __init__(self,graphsEdgeList, graphsNodeList):
        if len(graphsEdgeList) != len(graphsNodeList):
            raise Exception('the number of list in graphsEdgeList and graphsNodeList should be the same as the number of graphs')
        graphNumb = len(graphsEdgeList)
        self.graphs = []
        mergeEdges = []
        mergeNodes = []
        for i in range(graphNumb):   
            edges = graphsEdgeList[i]
            nodes = graphsNodeList[i]
            self.graphs.append(ValidGraph(edges, nodes))  
            mergeEdges.extend(edges)
            mergeNodes.extend(nodes)
      
        self.mergeGraph = ValidGraph(mergeEdges,mergeNodes) 
    
    def isConnect(self,upstream,downsteam,paramNum):
        if paramNum >=0 and paramNum < len(self.graphs):
            return self.graphs[paramNum].isConnect(upstream,downsteam)
        else:
            raise Exception("paramNum",paramNum,'is not available, max paramNum is',len(self.graphs)-1)


class GraphAll(GraphMng):
    def __init__(self,*args):
                #Smart-seq
        node1 = ['SRAToFastq',
                 'FastQC',
                 'FastqDump',
                 'Hisat2',
                 'SamToBam',
                 'Bamsort',
                 'BamSort',
                 'Stringtie',
                 'Tophat2',
                 'Star',
                 'Cufflinks',
                 'Cuffmerge',
                 'Cuffquant',
                 'Cuffnorm',
                 'Cuffdiff',
                 'HTSeq_sam2count',
                # 10x涉及的类的名称
                 'Quantification10x',
                 'PCA',
                # scATAC-seq
                 'FastQC',
                 'AdapterRemoval',
                 'Bowtie2',
                 'SamToBam',
                 'BamToSam',
                 'DuplicateRemoval',
                #Monocle
                 'MonocleQC',
                 'MonocleDC',
		 'Monocle2QC',
		 'Monocle2Pseudo',
                 #matrix preprocess
                 'matrixRw',     
                 # SC3
                 'SingleCellExperiment',
                 'SC3_DE',
                 'SC3_Cluster'            
                ]
        
        edge1 = [
                #Smart-seq
                ['FastqDump','Hisat2'],
                ['FastqDump','Tophat2'],
                ['FastqDump','FastQC'],
                ['Hisat2','HTSeq_sam2count'],
                ['Tophat2','HTSeq_sam2count'],
                ['Hisat2','SamToBam'],
                ['Tophat2','SamToBam'],
                ['SamToBam','Bamsort'],
                ['Bamsort','Cufflinks'],
                ['BamSort','Stringtie'],
                ['SamToBam','BamSort'],
                ['BamSort','Cufflinks'],
                ['Cufflinks','Cuffmerge'],
                ['Stringtie','Cuffmerge'],
                ['BamSort','Cuffquant'],
                ['Cuffquant','Cuffnorm'],
                ['Cuffmerge','Cuffdiff'],
                
            
                ['SRAToFastq','FastQC'],
                ['SRAToFastq','Tophat'],
                ['SRAToFastq','Star'],
                ['Tophat','Cufflinks'],
                ['Star','HTSeq'],
                #10x Genomeics
                ['Cellranger','Seuratpreprocessing'],
                ['Seuratpreprocessing','Seuratrun'],
                ['Qualification10x','PCA'],
                #drop-seq
                ['FastqToBam','BamMerge'],
                ['BamMerge','TagBarcode'],
                ['TagBarcode','TagBarcode'],
                ['TagBarcode','FilterBam'],
                ['FilterBam', 'TrimAdapter'],
                ['TrimAdapter','TrimPolyA'],
                ['TrimPolyA','BamToFastq'],
                ['BamToFastq','StarAlign'],
                ['StarAlign','SortBam'],
                ['TrimPolyA','MergeBamAlign'],
                ['MergeBamAlign','TagGene'],
                ['TagGene','DetectError'],
                ['DetectError','DigitalExpression'],
                ['DigitalExpression','MonocleQC'],
                ['DigitalExpression', 'EasyTreat'],
                ['MonocleQC','MonocleDC'],
                # ATAC-seq
                ['SRAToFastq', 'AdapterRemoval'],
                ['AdapterRemoval', 'Bowtie2'],
                ['Bowtie2', 'SamToBam'],
                ['SamToBam', 'BamSort'],
                ['BamSort', 'RmDuplicates'],
                ['RmDuplicates', 'LibComplexity'],
                ['RmDuplicates', 'BamToBed'],
                ['BamToBed', 'MergeToFrag'],
                ['BamToBed', 'RmChrOrMergeAllSample'],
                ['RmChrOrMergeAllSample', 'MergeToFrag'],
                ['MergeToFrag', 'FragLenDistri'],
                ['MergeToFrag', 'FragInPeak'],
                ['RmChrOrMergeAllSample', 'BedSort'],
                ['BedSort', 'PeakCalling'],
                ['PeakCalling', 'GenPeakWithFilter'],
                ['LibComplexity', 'CellFilter'],
                ['CellFilter', 'CellExtracterBam'],
                ['RmDuplicates', 'VarAndClustering'],
                ['CellExtracterBam', 'VarAndClustering'],

                #SC3
                ['SingleCellExperiment', 'SC3_DE'],
                ['SingleCellExperiment', 'SC3_Cluster'],
		#Monocle
		['Monocle2QC','Monocle2Pseudo'],
                ]

        edge2 = [
                ['SortBam', 'MergeBamAlign'],
                ['GenPeakWithFilter', 'VarAndClustering'],
                ['GenPeakWithFilter', 'FragInPeak'],
                ['FragInPeak', 'CellFilter'],
                ['RmDuplicates', 'CellExtracterBam'],
                ['Cuffmerge','Cuffquant'],
                ['Cuffmerge','Cuffnorm'],
                ['Cuffquant','Cuffdiff']
                ]

        super(GraphAll, self).__init__([edge1,edge2],[node1,[]])        
        
class GraphATACgl(GraphMng):
    def __init__(self,*args):
        edge1 = [
                ['UnzipAndMerge','AdapterRemoval'],
                ['UnzipAndMerge','FastQC'],
                ['AdapterRemoval','Bowtie2'],
                ['Bowtie2','SamToBam'],
                ['BamToSam','SamToBam'],
                ['SamToBam','BamToSam'],
                ]
        super(GraphATACgl, self).__init__([edge1],[[]])
