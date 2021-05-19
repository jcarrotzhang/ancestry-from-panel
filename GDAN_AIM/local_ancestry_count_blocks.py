#!/usr/bin/python -u
from __future__ import division
import sys, getopt
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

samplelist = str(sys.argv[1])
glban=open("gdan_aim_patient_ancestry_calls.txt")

eur_pop={}
afr_pop={}

for line in (raw.strip().split() for raw in glban):
       if "patient" not in line[0]:
               if "NA" not in line[len(line)-4]:
                       if (float(line[len(line)-4]) > 0.2) and (float(line[len(line)-4]) < 0.8) :
                              if (float(line[len(line)-5]) > 0.2) and (float(line[len(line)-5]) < 0.8):
                               if ("eur" in line[2]) or ("afr" in line[2]):
                                       if ("asian" not in line[2]) and ("sas" not in line[2]):
                                               (sampleID, tumor_type, consensus_call) = line[0:3]
                                               eur_pop["TCGA-"+str(sampleID)]=line[len(line)-5]
                                               afr_pop["TCGA-"+str(sampleID)]=line[len(line)-4]

bks=[]
OUTPUT=open("admixture_blocks.txt", 'w')

files=[]
for line in (raw.strip().split() for raw in samplelist):
	files.append(line[0])

for sample in eur_pop:
	print "working on", sample 
        inputfile1 = sample+"_A.bed"
        if inputfile1 in files:
        	input1=open(inputfile1) 
                for line1 in (raw.strip().split() for raw in input1):
                           if "UNK" not in line1[3]:
                                if str(line1[0])+":"+str(line1[1]) not in bks:
                                        bks.append(str(line1[0])+":"+str(line1[1]))
                                if str(line1[0])+":"+str(line1[2]) not in bks:
                                        bks.append(str(line1[0])+":"+str(line1[2]))
         inputfile2 = sample+"_B.bed"
         if inputfile2 in files:
                 input2=open(inputfile2)
                 for line1 in (raw.strip().split() for raw in input2):
                           if "UNK" not in line1[3]:
				if str(line1[0])+":"+str(line1[1]) not in bks:
                                        bks.append(str(line1[0])+":"+str(line1[1]))
                                if str(line1[0])+":"+str(line1[2]) not in bks:
					bks.append(str(line1[0])+":"+str(line1[2]))
                
for i in range(0, len(bks)):
	print >>OUTPUT, bks[i]
