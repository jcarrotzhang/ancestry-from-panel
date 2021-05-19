#!/usr/bin/python -u
from __future__ import division
import sys, getopt
import re
import resource
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
import statsmodels.formula.api as sm
import patsy
import statsmodels

inputfile = str(sys.argv[1])
loc = str(sys.argv[3])
Gene = str(sys.argv[2])
globan=open("gdan_aim_patient_ancestry_calls.txt")

cantype=("ACC","ESCA","LGG","PCPG","THCA","BLCA","GBM","LIHC","PRAD","THYM","BRCA","HNSC","LUAD","UCEC","CESC","KICH","LUSC","SARC","UCS","CHOL","KIRC","MESO","SKCM","UVM","COAD","KIRP","OV","STAD","DLBC","LAML","PAAD","TGCT")
chrom, pos = loc.split(":")
pos1, pos2 = pos.split("-")

gl=[]
eur_pop={}
afr_pop={}

### find admixed information. ###
for line in (raw.strip().split() for raw in globan):
	if "patient" not in line[0]:
		if "NA" not in line[len(line)-4]:
			if (float(line[len(line)-4]) > 0.1):
			#	if (float(line[len(line)-5]) < 0.9):
				if ("eur" in line[2]) or ("afr" in line[2]):
					if ("asian" not in line[2]) and ("sas" not in line[2]):
						(sampleID, tumor_type, consensus_call) = line[0:3]
						eur_pop["TCGA-"+str(sampleID)]=line[len(line)-5]
						afr_pop["TCGA-"+str(sampleID)]=line[len(line)-4]

### find local ancestry breakpoints and store in bks. ###
bks=[pos1, pos2]

samplelist=open(inputfile)
files=[]
for line in (raw.strip().split() for raw in samplelist):
               files.append(line[0])
for sample in eur_pop:
		inputfile1 = sample+"_A.bed"
		if inputfile1 in files:
                               input1=open(inputfile1) 
			       for line1 in (raw.strip().split() for raw in input1):
			      	       if (str(line1[0]) == chrom) and ("UNK" not in line1[3]):
					       if (int(line1[1]) < int(pos2)) and (int(line1[1]) > int(pos1)): 
						       if str(line1[1]) not in bks:
						    	   bks.append(str(line1[1]))
					       if (int(line1[2]) < int(pos2)) and (int(line1[2]) > int(pos1)): 
							if str(line1[2]) not in bks:
								bks.append(str(line1[2]))
		inputfile2 = sample+"_B.bed"
		if inputfile2 in files:
				input2=open(inputfile2)
				for line1 in (raw.strip().split() for raw in input2):
					if (str(line1[0]) == chrom) and ("UNK" not in line1[3]):
						if (int(line1[1]) < int(pos2)) and (int(line1[1]) > int(pos1)): 
							if str(line1[1]) not in bks:
								 bks.append(str(line1[1]))
						if (int(line1[2]) < int(pos2)) and (int(line1[2]) > int(pos1)):
							if str(line1[2]) not in bks:
								bks.append(str(line1[2]))

### find mutation status. ###								
mutation={}
for ctype in cantype:
	MC3 = "MC3_"+ctype+"_sample_sig_gene_table.txt"
	mc3=open(MCS)
	sample_pos={}
	first = mc3.readline().split('\t')
	for i in first:
		sampleID = i[:12]
		s= first.index(i)
		sample_pos[s]=sampleID
	for line in (raw.strip().split('\t') for raw in mc3):
		s = 1
		gene=str(line[0])
		if Gene in gene:
			for i in line[1:len(line)-1]:
				sid=sample_pos[s]
				if float(i) > 0:
					mutation[sid]=1
				else:
					mutation[sid]=0
				s = s + 1


### look for local ancestry for each bks. ####
lol={}

samplelist=open(inputfile)
files=[]
for line in (raw.strip().split() for raw in samplelist):
	files.append(line[0])
	
for sample in eur_pop:
			inputfile1 = sample+"_A.bed"
			if inputfile1 in files:
				input1=open(inputfile1)	
				for line2 in (raw.strip().split() for raw in input1):
					for key in bks:
						if key not in lol:
							lol[key]={}
						pos1=key
						if (str(line2[0]) == chrom) and (int(line2[1]) < int(pos1)) and (int(line2[2]) > int(pos1)):
							if "AFR" in line2[3]:
								a1=0
							if "EUR" in line2[3]:
								a1=1
							inputfile2 = sample+"_B.bed"
							if inputfile2 in files:
								input2=open(inputfile2)
								for line3 in (raw.strip().split() for raw in input2):
									if (str(line3[0]) == chrom) and (int(line3[1]) < int(pos1)) and (int(line3[2]) > int(pos1)):
										if "AFR" in line3[3]:
											if key in lol:
												if a1 == 0:
													lol[key][sample]=1
												if a1 == 1:
													lol[key][sample]=2
										if "EUR" in line3[3]:
											if key in lol:
												if a1 == 0:
													lol[key][sample]=0
												if a1 == 1:
													lol[key][sample]=1	
### perform association test. ####
for pos in lol:
	mu=[]
	an=[]
	gl_1=[]
	gl_2=[]
	for sample in lol[pos]:
		if sample in mutation:
			mu.append(int(mutation[sample]))
			an.append(int(lol[pos][sample]))
			gl_1.append(float(eur_pop[sample]))
			gl_2.append(float(afr_pop[sample]))

	df = pd.DataFrame({"mutation": mu, "ancestry": an, "eur":gl_1, "afr":gl_2})

	f = "mutation ~ ancestry + gl_1 + gl_2"
	y, X = patsy.dmatrices(f, df, return_type='dataframe')
	result = sm.Logit(y, X).fit()
	
	print pos, result.summary()

	f = "mutation ~ ancestry"
	y, X = patsy.dmatrices(f, df, return_type='dataframe')
	result = sm.Logit(y, X).fit()

	print pos, result.summary()

