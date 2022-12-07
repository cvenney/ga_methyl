#!/usr/bin/env Rscript
#13_plot_meth_vs_freqC.R

library(dplyr)

d <- read.table("10_DML_SNP_overlap_allelefreqs/allCpGs_cliff_pct_meth_pct_meth_freqC_DMLs_dwarf.txt", header=T)
n <- read.table("10_DML_SNP_overlap_allelefreqs/allCpGs_cliff_pct_meth_pct_meth_freqC_DMLs_normal.txt", header=T)

full <- inner_join(d, n)

write.table(full, "10_DML_SNP_overlap_allelefreqs/cliff_allCpGs_pct_meth_freqC_for_Eric.txt", row.names=FALSE, quote=FALSE)