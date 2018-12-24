# polygenic_selection

(A)
Use the polygenic_scores_pipeline to generate and plot polygenic scores for height using GIANT and UK Biobank GWAS for four 1000 Genomes European populations

# Update input parameters in the wrapper.sh script
Create and enter paths for the results and analysis folders
Enter the paths and names for the GWAS summary statistics files
Enter the paths and names for the VCF file (for example, 1000 genomes Phase 3) and the population label file (see 1000g_eur.txt as an example)

# Run as sh wrapper.sh
This script calls wrapper_present_day_samples.sh as well as other scripts for analysis and plotting

# Notes:
1) The input GWAS file should already be clumped (see plink --clump) or pruned for LD in a different way.
	# The provided GIANT and UKB files have been clumped using plink (r2 < 0.1, 1Mb, P < 0.01,  or r2 < 0.1, 1Mb, P < 5E-8)
	# Other GWAS summary statistics files can be used, simply update the paths to the GWAS files (GWASF1, GWASF2) in wrapper.sh. 
	  (Make sure the summary statistics files used follow the same column order as the example UKB and GIANT files provided) 
2) Only individuals (labelled by population) specified in POPFILE will be analyzed (see 1000g_eur.txt for an example)
3) If you are using a different prediction dataset than 1000 genomes Phase 3, 
	# update the vcfname variable in wrapper.sh to reflect the name of the vcf file
	# update the nums vector in compute_prs_indsubset_subpops_dg.py to reflect the sizes of the different populations 
	# update the x vector in visualize_prs_betase.r to reflect the names of the populations (as provided in POPFILE)
	# update the orderp and orderx vectors in plot_prs.r if using more than 4 populations 
4) If you have access to a cluster and multiple threads, you can run the loop in wrapper.sh in parallel on 22 threads to speed-up processing time
5) R library ggplot2 is required for the plotting scripts in the pipeline

(B)
Plot polygenic score figures (1000 Genomes and ancient DNA) from the Sohail, Maier et al, eLife (2019) manuscript

# Each figure folder provides the relevant plotting data and plotting script
Place the script and plotting data in the same directory and enter the path to the working directory in the mpath variable in the plotting script
Run the R script to create the plot, the plot will be output in the same folder 

# Notes:
1) Plotting scripts require ggplot2 library to be installed 
