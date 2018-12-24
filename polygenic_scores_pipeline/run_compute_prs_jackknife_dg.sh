#!/bin/bash 
### write out polygenic scores for different populations


VCFPATH="$1"
PREFIX="$2"
POPFILE="$3"
GWASP="$4"
i="$5"

GWASF="$PREFIX".effectacounts.indsubset.gwas
GWASFSP="$PREFIX".effectacounts.indsubset.subpops.gwas

echo -e "Step 2a: Combining all chromosomes"
cat $GWASP$PREFIX.*.effectacounts.indsubset.subpops.gwas.pvalthres.$i > $GWASP$PREFIX.effectacounts.indsubset.subpops.gwas.pvalthres.$i
	
echo -e "Step 2b: Computing polygenic score for each population"
python compute_prs_indsubset_subpops_dg.py $PREFIX $GWASP $GWASFSP.pvalthres.$i









