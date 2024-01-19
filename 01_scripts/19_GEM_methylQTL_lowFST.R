#!/usr/bin/env Rscript
#19_GEM_methylQTL.R
# conda activate R42
# srun -c 1 --mem 200G --time 7-00:00:00 -p medium -o 18_GEM_methylQTL_%j.log Rscript ./01_scripts/19_GEM_methylQTL_lowFST.R &

# methylQTL analysis -trying --mem 40G -c 8

# first three need columns (samples) in the same order!
meth="11_methylQTL/cliff_pct_meth_matrix_dml_freq.txt"
snp="11_methylQTL/all_maf0.05_pctind0.75_maxdepth25_cutOff0.7.geno.50pctlowestFST.GEM"
meta="11_methylQTL/cliff_covariates.txt"
pval=5.0E-08
output="11_methylQTL/cliff_methylQTL_results_50pctlowestFST.txt"

library(GEM)

#GEM_Gmodel(snp_file_name, covariate_file_name, methylation_file_name, Gmodel_pv, Gmodel_result_file_name)
GEM_Gmodel(snp, meta, meth, pval, output)

