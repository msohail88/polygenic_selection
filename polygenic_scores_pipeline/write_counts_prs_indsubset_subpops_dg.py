#!/usr/bin/env python

### write_counts_prs_indsubset_subpops_dg.py

import sys
import gzip
import random


mypath = sys.argv[1]
prefixf = sys.argv[2]
chrnum = sys.argv[3]
popfile = sys.argv[4]
vcffilename = sys.argv[5]
outfilename = prefixf + "." + chrnum + ".altcounts.indsubset.subpops"

vcf = gzip.open(mypath + vcffilename, "r")
outf = open(mypath + outfilename, "w")

### create individual dictionary #####
inddict = {}
with open(popfile, "r") as indf:
	next(indf)  #skip header
	for line in indf:
		mytabs = line.split()
		ind = mytabs[0]
		assig = mytabs[3]
		inddict[ind] = assig
		
subpops = set(inddict.values())
nsubpops = len(subpops)
subpopsl = list(subpops)
#print subpopsl



### open file and assess each snp ####
for line in vcf:
	if line.startswith("##"):
		continue
	if line.startswith("#CHROM"):
		headers = line.split()
		ids = headers[9:]
		totinds = len(ids)
		
		## check if in individual subset ###
		idsubset = []
		for id in ids:
			if id in inddict:
				idsubset.append(id)
		
		outheaderline = "\t".join(idsubset) + "\n"
		outheaderline2 = "\t".join(subpopsl) + "\n"

		outf.writelines(outheaderline)
		outf.writelines(outheaderline2)
		continue
	tabs = line.split()
	gts = tabs[9:]
	chrn = tabs[0]
	pos = tabs[1]
	rsid = tabs[2]
	refa = tabs[3]
	alta = tabs[4]
	
	# break totsnps into subpopulations
	tot = [0] * nsubpops
	alt = [0] *nsubpops
	
	outline = chrn + "\t" + pos+ "\t" + rsid + "\t" + refa + "\t" + alta
	### count number of alternative alleles and total number of alleles (given some alleles have missing data, some are haploid and some diploid) #######
	for i in range(0, totinds):
		gt = gts[i]
		id = ids[i]
		### only process if individual in id subset ########
		if id not in inddict:
			continue
		else:
			label = inddict[id]
			labeli = subpopsl.index(label)
			#print label
			#print str(labeli)
			if gt == "0|0":
				alt[labeli] += 0
				tot[labeli] += 2
			elif gt == "0|1" or gt == "1|0":
				alt[labeli] += 1
				tot[labeli] += 2
			elif  gt == "1|1":
				alt[labeli] += 2
				tot[labeli] += 2
			else:  ## missing
				pass

	### write out count of alt alleles and total alleles ####
	altout = ""
	for i in alt:
		altout = altout + "\t" + str(i)
	outline = outline + altout
		
	totout = ""
	for j in tot:
		totout = totout  + "\t" + str(j)
	outline = outline + totout + "\n"
	outf.writelines(outline)
	
vcf.close()
outf.close()




