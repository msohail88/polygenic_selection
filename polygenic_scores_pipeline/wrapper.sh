#!/bin/sh

##### INPUT PARAMETERS #####
# results folder
OUTPATH="./results/"
# P-value threshold for polygenic score computation 
pval=0.01
## gwas summary statistics 1
# gwas summary statistics file
GWASF1=GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.header.txt.clumpedout.0.01.gz 
# prefix for output files
GWASPREFIX1="GIANT"
# output folder for intermediate files from first gwas analysis
GWASP1="./gwas1/analyses/"
## gwas summary statistics 2
# gwas summary statistics file
GWASF2=50.assoc.tsv.processed.nodups.clumpedout.0.01.gz 
# prefix for output files
GWASPREFIX2="UKB"
# output folder for intermediate files from second gwas analysis
GWASP2="./gwas2/analyses/"
## vcf file and population labels file
VCFPATH="./1000Genomes_Phase3/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
# population labels file
POPFILE="$VCFPATH"1000g_eur.txt
##### End of INPUT PARAMETERS #####

## write out count files for alternate alleles for a given list of individuals outputting counts for each population
echo "PART A: Runnning script to obtain allelic counts from vcf files for all populations specified (this is the most time-consuming part)"
PREFIX="ALL"
#loopStart,i    
for i in `seq 22`; do            
    printf "Processing chr %s\n" $i
    vcfname="$PREFIX".chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    python write_counts_prs_indsubset_subpops_dg.py $VCFPATH $PREFIX "$i" $POPFILE $vcfname
#loopEnd                     
done

# run wrapper script to determine which allele is the risk alleles for GWAS SNPs and compute polygenic scores using gwas summary statistics
echo "PART B: Runnning scripts to determine the risk allele for GWAS SNPs and compute polygenic scores using different GWAS summary statistics"
sh wrapper_present_day_samples.sh $GWASF1 $GWASPREFIX1 $pval $GWASP1 $OUTPATH $VCFPATH $POPFILE
sh wrapper_present_day_samples.sh $GWASF2 $GWASPREFIX2 $pval $GWASP2 $OUTPATH $VCFPATH $POPFILE

# combine the two GWAS results into a single file for plotting
echo "PART C: Preparing output data frame for plotting"
Rscript make_output_table_betase.r $OUTPATH $pval $GWASPREFIX1 $GWASPREFIX2

# Make plot 
echo "PART D: Creating plot showing polygenic scores using different GWAS"
Rscript plot_prs.r $OUTPATH $pval

