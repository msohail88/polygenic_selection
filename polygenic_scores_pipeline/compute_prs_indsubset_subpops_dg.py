#!/usr/bin/env python

### compute_prs_indsubset_subpops_dg.py

import sys
import gzip


prefixf = sys.argv[1]
gwaspath = sys.argv[2]
gwasfile = sys.argv[3]
outfilename = gwasfile + ".prs"
outfilename2 = gwasfile + ".prs.noheader"

outf = open(gwaspath + outfilename, "w")
outf2 = open(gwaspath + outfilename2, "w")

nums = [99,107,91,107]  #'CEU', 'IBS', 'GBR', 'TSI'


tot = sum(nums)
weights = [0] * len(nums)
for i in range(0, len(nums)):
	weights[i] = nums[i]/float(tot)


# read in counts file and write out prs
with open(gwaspath + gwasfile, "r") as f:
	header1 = next(f)
	header2 = next(f)
	outf.writelines(header1)
	outf.writelines(header2)
	
	# compute prs
	nsubpops = len(header2.split("\t"))
	prs = [0] * nsubpops
	tots = [0] * nsubpops
	counts = [0] * nsubpops
	varprs = [0] * nsubpops
	mean = 0.0 
	Va = 0.0
	for line in f:
		if line == header1 or line == header2:
			continue 
		tabs = line.split("\t")
		flip = int(tabs[(nsubpops*2) + 5])
		b = float(tabs[(nsubpops*2) + 6])
		avgpsum = 0.0
		for pi in range(0, nsubpops):
			ci = float(tabs[5+pi])
			ti = float(tabs[5+nsubpops+pi])
			if ti > 0: 
				prs[pi] += (ci/ti) * b
				tots[pi] += ti
				counts[pi] += ci
				avgpsum = avgpsum + ((ci/ti)*weights[pi])
				## for beta standard errors
				alpha = 1 + ci
				beta = 1 + ti - ci
				num = alpha * beta
				denom = ((alpha+beta)**2)*(alpha+beta+1)
				varpl = num /denom
				varprs[pi] += (b**2) * varpl
		avgp = avgpsum
		mean = mean + (b * avgp)
		#print str(b), str(avgp)
		Va = Va + ((b**2) * avgp * (1 - avgp))
	
	
			
	# create outline
	#avgs = [0] * nsubpops 
	outline = ""
	
	# compute overall mean and VA
	meanout = mean
	Vaout = Va
	
	# add sums to outline, compute avgs
	for pi in range(0, nsubpops):
		outline += str(prs[pi]) + "\t"
		#avgs[pi] = prs[pi]/counts[pi]
	
	# add tot non-missing alleles to line
	for pi in range(0, nsubpops):
		outline += str(tots[pi]) + "\t"
	
	# add tot avgs to line 
	for pi in range(0, nsubpops):
		outline += str(counts[pi]) + "\t"
	
	
	outline = outline + str(meanout) + "\t" + str(Vaout) + "\t"
	
	# add beta variance of polygenic score for each population
	for pi in range(0, nsubpops):
		outline += str(varprs[pi]) + "\t"	
	
	outline = outline + "\n"
	outf.writelines(outline)
	outf2.writelines(outline)
	
outf.close()
outf2.close()
f.close()
		
		
		
		


