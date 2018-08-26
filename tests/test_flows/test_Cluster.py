#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2018-04-07 10:27:03
# @Author  : Zhenyi Wang
# @Email: wangzy17@mails.tsinghua.edu.cn

from hcacn.core import Configure,Schedule
from hcacn.flows import Cluster

Configure.setIdentity('Clustertest')
#test monocle 
clusterObj= Cluster(matrixdata="./minidata/Monocle/out_gene_exon_tagged.dge.txt",
	                algorithm = 'Monocle', resultDir='./result')()
#test sc3
#ClusterObj= Cluster(sc3matrix_file = "./minidata/sc3cluster/",
#	                sc3ann = './minidata/sc3annotation/',
#	                algorithm = 'sc3', resultDir='./result')()