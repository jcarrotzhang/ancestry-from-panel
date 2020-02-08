#!/usr/bin/python -u
from __future__ import division
#import sys, getopt
import re
import numpy as np
import argparse, os, sys
from scipy import stats
import math
import matplotlib
import matplotlib.pyplot as plt
matplotlib.pyplot.switch_backend('agg')
from scipy.stats import mannwhitneyu
import pandas as pd
import patsy

cantype=("ACC","ESCA","LGG","PCPG","THCA","BLCA","GBM","LIHC","PRAD","THYM","BRCA","HNSC","LUAD","UCEC","CESC","KICH","LUSC","SARC","UCS","CHOL","KIRC","MESO","SKCM","UVM","COAD","KIRP","OV","STAD","DLBC","LAML","PAAD","TGCT")

input2=open("TCGA_geno/gdan_aim_patient_ancestry_calls.txt")
input3=open("/xchip/gcc_data/results2/aggregate/sif/production.postprocess.sif.txt")

eur_pop={}
afr_pop={}

arrayids={}
for line in (raw.strip().split() for raw in input3):
       m = re.search('GenomeWideSNP_6_(.+?)$', line[3])
       if m:
       	arrayid=m.group(1)
	samplename=str(line[2])
	arrayids[samplename]=arrayid

for line in (raw.strip().split() for raw in input2):
       if "patient" not in line[0]:
               if "NA" not in line[len(line)-4]:
                       if (float(line[len(line)-4]) > 0.2) and (float(line[len(line)-4]) < 0.8) :
                              if (float(line[len(line)-5]) > 0.2) and (float(line[len(line)-5]) < 0.8):
                               if ("eur" in line[2]) or ("afr" in line[2]):
                                       if ("asian" not in line[2]) and ("sas" not in line[2]):
                                               (sampleID, tumor_type, consensus_call) = line[0:3]
                                               eur_pop["TCGA-"+str(sampleID)]=line[len(line)-5]
                                               afr_pop["TCGA-"+str(sampleID)]=line[len(line)-4]
print len(eur_pop)

bks=[]
OUTPUT=open("TCGA_geno/analysis/admixture_blocks.txt", 'w')
for ctype in cantype:
       samplelist=open("/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/samplelist")
       files=[]
       for line in (raw.strip().split() for raw in samplelist):
               files.append(line[0])
       for sample in eur_pop:
               if sample in arrayids:
                       sn=arrayids[sample]
		       print "working on", ctype, sample 
                       inputfile4 = "/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/"+sn+"_A.bed"
                       if inputfile4 in files:
                              input4=open(inputfile4) 
                              for line1 in (raw.strip().split() for raw in input4):
                                      if "UNK" not in line1[3]:
                                                      if str(line1[0])+":"+str(line1[1]) not in bks:
                                                          bks.append(str(line1[0])+":"+str(line1[1]))
                                                      if str(line1[0])+":"+str(line1[2]) not in bks:
                                                          bks.append(str(line1[0])+":"+str(line1[2]))
                       inputfile5 = "/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/"+sn+"_B.bed"
                       if inputfile5 in files:
                               input5=open(inputfile5)
                               for line1 in (raw.strip().split() for raw in input5):
                                       if "UNK" not in line1[3]:
						       if str(line1[0])+":"+str(line1[1]) not in bks:
                                                                bks.append(str(line1[0])+":"+str(line1[1]))
                                                       if str(line1[0])+":"+str(line1[2]) not in bks:
								bks.append(str(line1[0])+":"+str(line1[2]))
                
for i in range(0, len(bks)):
	print >>OUTPUT, bks[i]
