#!/bin/sh

##### INPUT PARAMETERS #####
## Specify params for GWAS file
# Specify path and name of GWAS summary statistics file
GWASF="$1"
# Specify GWAS name for outputted analysis files, for example, "GIANT"
GWASPREFIX="$2"
# Specify P-value threshold
pval="$3"
# Specify output path for intermediate analysis files
GWASP="$4"
## Specify path for output plots
OUTPATH="$5"
## specify location of vcffile 
VCFPATH="$6"
## specify location and name of population file
POPFILE="$7"
##### End of INPUT PARAMETERS #####

printf "GWAS being processed: %s\n" $GWASPREFIX
## switch counts files for alternate alleles generated in the last part to count files for effect alleles in all given populations
printf "Step 1: Flipping alternate allele for effect allele for all GWAS SNPs at a given P-value threshold in %s\n" $GWASPREFIX
PREFIX="ALL"
#loopStart,i    
for i in `seq 22`; do            
    
    printf "Processing chr %s\n" $i
	python write_counts_prs_indsubset_subpops_riska_pruned_dg_forclumped.py $VCFPATH $PREFIX "$i" $GWASP $GWASF $pval

#loopEnd                     
done
    
# compute polygenic scores
printf "Step 2: Computing polygenic scores using %s\n" $GWASPREFIX
sh run_compute_prs_jackknife_dg.sh $VCFPATH $PREFIX $POPFILE $GWASP $pval

# write out plotting summaries of polygenic scores   
printf "Step 3: Outputting polygenic scores with confidence intervals for all populations using %s\n" $GWASPREFIX
Rscript visualize_prs_betase.r $GWASP $OUTPATH $POPFILE $PREFIX $pval $GWASPREFIX









