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

cantype=["ACC","ESCA","LGG","PCPG","THCA","BLCA","GBM","LIHC","PRAD","THYM","BRCA","HNSC","LUAD","UCEC","CESC","KICH","LUSC","SARC","UCS","CHOL","KIRC","MESO","SKCM","UVM","COAD","KIRP","OV","STAD","DLBC","LAML","PAAD","TGCT"]
loc = str(sys.argv[2])
Gene = str(sys.argv[1])
chrom, pos = loc.split(":")
pos1, pos2 = pos.split("-")

gl=[]

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

### find admixed information. ###
for line in (raw.strip().split() for raw in input2):
	if "patient" not in line[0]:
		if "NA" not in line[len(line)-4]:
			if (float(line[len(line)-4]) > 0.1):
			#	if (float(line[len(line)-5]) < 0.9):
				if ("eur" in line[2]) or ("afr" in line[2]):
					if ("asian" not in line[2]) and ("sas" not in line[2]):
						(sampleID, tumor_type, consensus_call) = line[0:3]
						eur_pop["TCGA-"+str(sampleID)]=line[len(line)-5]
						afr_pop["TCGA-"+str(sampleID)]=line[len(line)-4]

bks=[pos1, pos2]

for ctype in cantype:
       samplelist=open("/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/samplelist")
       files=[]
       for line in (raw.strip().split() for raw in samplelist):
               files.append(line[0])
       for sample in eur_pop:
               if sample in arrayids:
			#print "working on", ctype, sample 
			sn=arrayids[sample]
			inputfile4 = "/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/"+sn+"_A.bed"
			if inputfile4 in files:
                               input4=open(inputfile4) 
			       for line1 in (raw.strip().split() for raw in input4):
			               if (str(line1[0]) == chrom) and ("UNK" not in line1[3]):
					       if (int(line1[1]) < int(pos2)) and (int(line1[1]) > int(pos1)): 
						       if str(line1[1]) not in bks:
						    	   bks.append(str(line1[1]))
					       if (int(line1[2]) < int(pos2)) and (int(line1[2]) > int(pos1)): 
							if str(line1[2]) not in bks:
								bks.append(str(line1[2]))
			inputfile5 = "/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/"+sn+"_B.bed"
			if inputfile5 in files:
				input5=open(inputfile5)
				for line1 in (raw.strip().split() for raw in input5):
					if (str(line1[0]) == chrom) and ("UNK" not in line1[3]):
						if (int(line1[1]) < int(pos2)) and (int(line1[1]) > int(pos1)): 
							if str(line1[1]) not in bks:
								 bks.append(str(line1[1]))
						if (int(line1[2]) < int(pos2)) and (int(line1[2]) > int(pos1)):
							if str(line1[2]) not in bks:
								bks.append(str(line1[2]))

mutation={}
for ctype in cantype:
	inputfile1 = "TCGA_geno/mc3"+"/MC3_"+ctype+"_sample_sig_gene_table.txt"
	input1=open(inputfile1)
	sample_pos={}
	first = input1.readline().split('\t')
	for i in first:
		sampleID = i[:12]
		s= first.index(i)
		sample_pos[s]=sampleID
	for line in (raw.strip().split('\t') for raw in input1):
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


### look for local ancestry. ####
lol={}
for ctype in cantype:
	samplelist=open("/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/samplelist")
	files=[]
	for line in (raw.strip().split() for raw in samplelist):
		files.append(line[0])
	for sample in eur_pop:
		if sample in arrayids:
			#print "working on", ctype, sample 
			sn=arrayids[sample]
			inputfile4 = "/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/"+sn+"_A.bed"
			if inputfile4 in files:
				input4=open(inputfile4)	
				for line2 in (raw.strip().split() for raw in input4):
					for key in bks:
						if key not in lol:
							lol[key]={}
						pos1=key
						if (str(line2[0]) == chrom) and (int(line2[1]) < int(pos1)) and (int(line2[2]) > int(pos1)):
							if "AFR" in line2[3]:
								a1=0
							if "EUR" in line2[3]:
								a1=1
							inputfile5 = "/cga/meyerson/home/zhangj/kgp/TCGA_geno/0205_seunghun/"+ctype+"/final/nofbk/"+sn+"_B.bed"
							if inputfile5 in files:
								input5=open(inputfile5)
								for line3 in (raw.strip().split() for raw in input5):
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

