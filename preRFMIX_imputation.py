#!/usr/bin/python -u
from __future__ import division
import sys, getopt
import re
import argparse, os, sys
import gzip

case = str(sys.argv[1])
ref = str(sys.argv[2])
loc = str(sys.argv[3])
alleles = str(sys.argv[4])
snp_location = str(sys.argv[5])+".snp_locations"
maps = str(sys.argv[5])+".map"
classes = str(sys.argv[6])
chrom = str(sys.argv[7])
pop=open("1000GP_Phase3.sample")
eur_s=[]
eas_s=[]
afr_s=[]

for line in (raw.strip().split(' ') for raw in pop):
        if ("EUR" in line[2]) and ("FIN" not in line[1]):
                eur_s.append(line[0])
        if ("EAS" in line[2]):
                eas_s.append(line[0])
        if ("AFR" in line[2]):
                afr_s.append(line[0])

eur_samples=eur_s[0:500]
eas_samples=eas_s[0:500]
afr_samples=afr_s[0:500]
n_eas=["1"]*(len(eas_samples)*2)
n_eur=["2"]*(len(eur_samples)*2)
n_afr=["3"]*(len(afr_samples)*2)

snp=open("random_SNPs.txt")
rids=[]
rr={}
for line in (raw.strip().split('\t') for raw in snp):
        if line[0] == chrom:
                rids.append(str(line[1]))
                rr[str(line[1])]=line[2]

location=open(loc)
ll={}
for line in (raw.strip().split(' ') for raw in location):
        if str(line[0]) in rids:
                ll[str(line[0])]=str(line[2])

pp={}
input2=open(case)
n = 0
for line in (raw.strip().split('\t') for raw in input2):
        if "#" not in line[0]:
                rrid = str(line[1])
                if rrid in ll:
                        if n == 2:
                                n = 0
                                phasing=line[9:len(line)]
                                for i in range(0, len(phasing)):
                                        l = phasing[i].split(":")[0]
                                        p="".join(l.split("|"))
                                        if ("10" in p) or ("01" in p) or ("00" in p) or ("11" in p): ### exclude biallelic SNPs. ###
                                                if rrid in pp:
                                                        pp[rrid] = pp[rrid]+""+str(p)
                                                else:
                                                        pp[rrid] = str(p)
                        else:
                                n = n + 1

        if "#CHROM" in line[0]:
                n_case=["0"] * ((len(line)-9) * 2)


OUT=open(classes, 'w')
out=" ".join(n_case)+" "+" ".join(n_eas)+" "+" ".join(n_eur)+" "+" ".join(n_afr)

print >>OUT, out
input1=gzip.open(ref)
pp_eur={}
pp_afr={}
pp_eas={}
for line in (raw.strip().split('\t') for raw in input1):
        if "#CHROM" in line[0]:
                sampleName = line[9:len(line)]
        if ("#" not in line[0]) and ("#CHROM" not in line[0]):
                       rid = str(line[1])
                       if rid in pp:
                                phasing=line[9:len(line)]
                                for i in range(0, len(sampleName)):
                                        l = phasing[i].split(":")[0]
                                        p="".join(l.split("|"))
                                        if ("10" in p) or ("01" in p) or ("00" in p) or ("11" in p): ### exclude biallelic SNPs. ###
                                                if sampleName[i] in eas_samples:
                                                        if rid in pp_eas:
                                                                pp_eas[rid] = pp_eas[rid]+""+str(p)
                                                        else:
                                                                pp_eas[rid] = str(p)
                                                if sampleName[i] in afr_samples:
                                                        if rid in pp_afr:
                                                                pp_afr[rid] = pp_afr[rid]+""+str(p)
                                                        else:
                                                                pp_afr[rid] = str(p)
                                                if sampleName[i] in eur_samples:
                                                        if rid in pp_eur:
                                                                pp_eur[rid] = pp_eur[rid]+""+str(p)
                                                        else:
                                                                pp_eur[rid] = str(p)


OUT1=open(alleles, 'w')
OUT2=open(snp_location, 'w')
OUT3=open(maps, 'w')
for p in sorted(pp.iterkeys(), key=int):
        if (p in pp_eas) and  (p in pp_eur) and  (p in pp_afr):
                out=pp[p]+""+pp_eas[p]+""+pp_eur[p]+""+pp_afr[p]
                #out=pp[p]+""+pp_eas[p]+""+pp_eur[p]
                out2=str(p)+" "+str(ll[p])+" "+str(rr[p])
                print >>OUT1, out
                print >>OUT2, ll[p]
                print >>OUT3, out2
