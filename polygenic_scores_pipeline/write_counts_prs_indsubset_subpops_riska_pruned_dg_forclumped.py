#!/usr/bin/env python

### write_counts_prs_indsubset_subpops_riska_pruned_dg_forclumped.py

import sys
import gzip

mypath = sys.argv[1]
prefixf = sys.argv[2]
chrnum = sys.argv[3]
gwaspath = sys.argv[4]
gwasfile = sys.argv[5]
pvalthres = sys.argv[6]

countsfilename = prefixf + "." + chrnum + ".altcounts.indsubset.subpops"
outfilename = prefixf + "." + chrnum + ".effectacounts.indsubset.subpops.gwas.pvalthres." + pvalthres

outf = open(gwaspath + outfilename, "w")

 

## store snp information for all GWAS snps
snpdict = {}
with gzip.open(gwasfile, "r") as gwasf:
	next(gwasf) ## skip header
	for line in gwasf:
		tabs = line.strip().split("\t")
		rsid = tabs[0]
		a1 = tabs[1]
		a2 = tabs[2]
		b = tabs[4]
		try:
			p = float(tabs[6])
		except:
			continue
		if p < float(pvalthres):
			snpdict[rsid] = [a1,a2,b]	
	
## read in alta counts file and write out effecta counts file
with open(mypath + countsfilename) as f:
	# skip 1st and 2nd header lines
        header1 = next(f) 
        header2 = next(f)
        outf.writelines(header1)
	outf.writelines(header2)
	subpops = header2.split()
	nsubpops = len(subpops)
	# go through all snps
	for line in f: 
		tabs = line.split("\t")
		chr = tabs[0]
		pos = tabs[1]
		rsid = tabs[2]
		refa = tabs[3]
		alta = tabs[4]
		counts = [0] * nsubpops
		tots = [0] * nsubpops
		counter = 0
		for ni in range(0, nsubpops):
			ci = int(tabs[5+ni])
			ti = int(tabs[5+nsubpops+ni])
			counts[counter] = ci
			tots[counter] = ti 
			counter = counter + 1 

		if rsid in snpdict:
			params = snpdict[rsid] 
			a1 = params[0]
			a2 = params[1]
			b = params[2]
			# check if effect allele is alt allele 
			flip = "0"
			if a1 == alta:
				outline = line.strip() + "\t" + flip + "\t" + b +  "\n"
			else:
				flip = "1"
				outc = [0] * nsubpops
				outline = chr + "\t" + pos + "\t" + rsid + "\t" + refa + "\t" + alta  + "\t" 
				for i in range(0, nsubpops):
					outci = tots[i] - counts[i]
					outc[i]=outci
					outline = outline + str(outci) + "\t" 
				for i in range(0, nsubpops):
					outline = outline + str(tots[i]) + "\t" 
				outline = outline + flip + "\t" + b +"\n"
			
			# write out outline
			outf.writelines(outline)

f.close()
gwasf.close()
			
		
		



