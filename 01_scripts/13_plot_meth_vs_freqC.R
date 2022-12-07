#!/usr/bin/env Rscript
#13_plot_meth_vs_freqC.R

library(dplyr)

d <- read.table("10_DML_SNP_overlap_allelefreqs/cliff_pct_meth_pct_meth_freqC_DMLs_dwarf.txt", header=T)
n <- read.table("10_DML_SNP_overlap_allelefreqs/cliff_pct_meth_pct_meth_freqC_DMLs_normal.txt", header=T)

full <- inner_join(d, n)

pdf(file = "10_DML_SNP_overlap_allelefreqs/cliff_dwarf_scatterplot.pdf") 
# plot
	plot(full$pct_meth_dwarf, full$freq_C_dwarf, main="freqC x meth Cliff",
		xlab="%M", ylab="freqC", pch=1, col="midnightblue") 
	points(full$pct_meth_normal, full$freq_C_normal,
		xlab="%M", ylab="freqC", pch=1, col="gold3")
dev.off()