#!/usr/bin/pythoni -u
from __future__ import division
import sys, getopt
import re
import pysam
#import pysamstats
	
base_quality = 20
mapping_quality = 20

reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
#t_samfile = pysam.AlignmentFile("Meyerson_LA_NSCLC_B5_SP30006e/P4765-p5-B7_L_160716-CEPH1328-p1-H1_L_159889/P4765-p5-B7_L_160716.dedup.cleaned.bam")
#n_samfile = pysam.AlignmentFile("Meyerson_LA_NSCLC_B5_SP30006e/P4765-p5-B7_L_160716-CEPH1328-p1-H1_L_159889/CEPH1328-p1-H1_L_159889.dedup.cleaned.bam")
tumorfile = sys.argv[1]
bed = "ANALYSIS/EGFR_KRAS_TP53_hotspots.bed"
output = open("ANALYSIS/EGFR_KRAS_TP53_scriptcalled_v3.txt", 'w')
ref = pysam.FastaFile(reference)

samplelist=open(tumorfile)
for line in (raw.strip().split() for raw in samplelist):
	sample = line[0]
	t_samfile = pysam.AlignmentFile(sample)
	
	bedfile=open(bed)
	array=[]
	for bedline in (raw.strip().split() for raw in bedfile):
		chrom=str(bedline[0])
		pos=int(bedline[1])

		if ("EGFR_del_19" in bedline[2]) or ("EGFR_ins" in bedline[2]) or ("EGFR_del_18" in bedline[2]) or ("MET_Splice" in bedline[2]) or ("ERBB2_ins" in bedline[2]):
			for pileupcolumn in t_samfile.pileup(chrom, pos-10, pos+10):
				if (pileupcolumn.pos < pos+10) and (pileupcolumn.pos > pos-10):
					alt = 0; ref = 0
					for pileupread in pileupcolumn.pileups:
						if pileupread.alignment.is_proper_pair and not pileupread.alignment.is_duplicate:
							if pileupread.alignment.mapping_quality >= mapping_quality:
								if pileupread.is_del or pileupread.is_refskip or pileupread.indel:
									alt = alt + 1
								else:
									ref = ref + 1

					if alt+ref > 5:
					#if (alt == 2) and (ref < 55):
					#	if sample not in array:
					#		out=str(bedline[2])+"\t"+str(sample)+"\t"+str(alt)+"\t"+str(ref)
					#		print out
					#		array.append(sample)

						if (alt > 2) and (alt/(alt+ref) > 0.05):
							if sample not in array:
								out=str(bedline[2])+"\t"+str(sample)+"\t"+str(alt)+"\t"+str(ref)
								print >>output,  out
								array.append(sample)
								break

		#if bedline[2] in ["EGFR_L858R", "EGFR_L861Q", "EGFR_S768I", "BRAF_V600E"]:
		if ("EGFR_L858R" in bedline[2]) or ("EGFR_L861Q" in bedline[2]) or ("EGFR_S768I" in bedline[2]) or ("BRAF_V600E" in bedline[2]):	
			alt=str(bedline[3])
			for pileupcolumn in t_samfile.pileup(chrom, pos-1, pos):
				if pileupcolumn.pos == pos-1:
					ALT = 0; ref = 0
					for pileupread in pileupcolumn.pileups:
						if pileupread.alignment.is_proper_pair and not pileupread.alignment.is_duplicate:
							if not pileupread.is_del and not pileupread.is_refskip:
								if pileupread.alignment.mapping_quality >= mapping_quality and pileupread.alignment.query_qualities[pileupread.query_position] >= base_quality:
									altseq = pileupread.alignment.query_sequence[pileupread.query_position]
									if altseq == alt:
										ALT += 1
									else:
										ref += 1
			
			if ALT+ref > 5:
				#	if (ALT == 2) and (ref < 20):
				#		if sample not in array:
				#			out=str(bedline[2])+"\t"+str(sample)+"\t"+str(ALT)+"\t"+str(ref)
				#			print out
				#			array.append(sample)
			
				if (ALT > 2) and (ALT/(ALT+ref) > 0.05):
					if sample not in array:
						out=str(bedline[2])+"\t"+str(sample)+"\t"+str(ALT)+"\t"+str(ref)
						print >>output,  out
						array.append(sample)

		#if bedline[2] in ["KRAS_G12N", "KRAS_G13N", "KRAS_Q61N", "KRAS_A146N", "BRAF_G466N", "BRAF_G469N", "MET_Splice", "EGFR_K754N", "EGFR_G719N"]:
		if ("KRAS_G12N" in bedline[2]) or ("BRAF_G469N" in bedline[2]) or ("EGFR_G719N" in bedline[2]) or ("KRAS_G13N" in bedline[2]) or ("KRAS_Q61N" in bedline[2]) or ("KRAS_A146N" in bedline[2]) or ("BRAF_G466N" in bedline[2]) or ("MET_Splice" in bedline[2]) or ("EGFR_K754N" in bedline[2]) or ("EGFR_G719N"in bedline[2])  or ("TP53" in bedline[2]):
			
			REF=str(bedline[3])
			for pileupcolumn in t_samfile.pileup(chrom, pos-1, pos):
				if pileupcolumn.pos == pos-1:
					ALT = 0; ref = 0
					for pileupread in pileupcolumn.pileups:
						if pileupread.alignment.is_proper_pair and not pileupread.alignment.is_duplicate:
							if not pileupread.is_del and not pileupread.is_refskip:
								if pileupread.alignment.mapping_quality >= mapping_quality and pileupread.alignment.query_qualities[pileupread.query_position] >= base_quality:
									altseq = pileupread.alignment.query_sequence[pileupread.query_position]
									if altseq is not REF:
										ALT += 1
									else:
										ref += 1

			if ALT+ref > 5:
					#if (ALT == 2) and (ref < 20):
					#	if sample not in array:
					#		out=str(bedline[2])+"\t"+str(sample)+"\t"+str(ALT)+"\t"+str(ref)
					#		print out
					#		array.append(sample)

				if (ALT > 2) and (ALT/(ALT+ref) > 0.05):
					if sample not in array:
						out=str(bedline[2])+"\t"+str(sample)+"\t"+str(ALT)+"\t"+str(ref)
						print >>output,  out
						array.append(sample)















