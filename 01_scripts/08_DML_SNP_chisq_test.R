#!/usr/bin/env Rscript
#08_DML_SNP_chisq_test.R
# conda activate R411 before use
#srun -c 1 --mem 80G -p small --time 1-00:00:00 -J 08_DML_SNP_test -o 08_DML_SNP_chisq_test_%j.log Rscript ./01_scripts/08_DML_SNP_chisq_test.R &

#use "gunzip -c 03_SNP_data/all_pctind0.5_maxdepth25.bed.gz | wc -l" for genome length
GEN_LENGTH=1998994058

library(regioneR)

DML <- toGRanges("08_DSS_results/cliff_DMLs_p0.05.bed")
CG <- toGRanges("/project/lbernatchez/users/clven7/whitefish/am/ga_permutation_cliff/04_permutation_files/all_CpG_covered_by_SNPs_nonamecheck_uniques_3col.bed")
outlierSNP <- toGRanges("02_SNPs/CN_CD_maf0.05_pctind0.75_maxdepth25_bypos_5pct.bed")
allSNP <- toGRanges("/project/lbernatchez/users/clven7/whitefish/am/ga_permutation_cliff/03_SNP_data/all_maf0.05_pctind0.75_maxdepth25.bed.gz")
CPGwSNPcov <- toGRanges("08_DSS_results/cliff_DMLs_all_CpGs.bed")


### chi-sq based on CpG sites only
print("####################################### chi-sq based on CpG sites only #######################################")
nCpG <- as.numeric(length(CG))
print(paste("number of CpG sites covered by SNP data =", nCpG))
print(paste("number of CpG site nucleotides covered by SNP data =", nCpG*2))

nDML <- as.numeric(length(DML))
print(paste("number of DMLs =", nDML))

nOutlier <- as.numeric(length(outlierSNP))
print(paste("number of outlier SNPs in analysis =", nOutlier))

overlaps <- as.numeric(numOverlaps(DML, outlierSNP))
print(paste("number of DML-SNP overlaps =", overlaps))

exp_overlap <- (nDML/nCpG) * (nOutlier/nCpG) * 2
print(paste("expected DML-SNP overlap rate", exp_overlap))
obs_overlap <- overlaps / nCpG
print(paste("observed DML-SNP overlap rate", obs_overlap))

nCpG2 <- (nCpG*nCpG)
nDML_nOutlier <- nDML*nOutlier*2

cont_tab<-rbind(c((nCpG2),(nDML_nOutlier)),
                  c(nCpG,overlaps))
print(paste("contingency table", cont_tab))

x2 <- chisq.test(cont_tab, correct=FALSE)
print(paste("chi-sq results:"))
print(x2)


### chi-sq based on entire genome
print("####################################### chi-sq based on entire genome #######################################")
exp_overlap_wg <- (nDML/GEN_LENGTH) * (nOutlier/GEN_LENGTH) * 2
print(paste("expected DML-SNP overlap rate", exp_overlap_wg))
obs_overlap_wg <- overlaps/GEN_LENGTH
print(paste("observed DML-SNP overlap rate", obs_overlap_wg))

GEN_LENGTH2 <- as.numeric(GEN_LENGTH * GEN_LENGTH)

cont_tab<-rbind(c((GEN_LENGTH2),(nDML_nOutlier)),
                  c(GEN_LENGTH,overlaps))
print(paste("contingency table", cont_tab))

x2 <- chisq.test(cont_tab, correct=FALSE)
print(paste("chi-sq results:"))
print(x2)


### chi-sq based on CpGs - only considering SNPs in CpGs
print("####################################### chi-sq based on only outlier SNPs in CpGs #######################################")
#nDML/nCpG * nOutliersInCpG/nCpG= ???
nCpGoutlier <- as.numeric(numOverlaps(CG, outlierSNP))
print(paste("number of outlier SNPs in CpG sites", nCpGoutlier))

exp_overlap_ocg <- (nDML/nCpG) * (nCpGoutlier/nCpG)
print(paste("expected DML-outlier CpG site SNP overlaps", exp_overlap_ocg))
obs_overlap <- overlaps / nCpG
print(paste("observed DML-SNP overlap rate", obs_overlap))

cont_tab<-rbind(c((nCpG2),((nCpGoutlier*nDML))),
                  c(nCpG,overlaps))
print(paste("contingency table", cont_tab))

x2 <- chisq.test(cont_tab, correct=FALSE)
print(paste("chi-sq results:"))
print(x2)


### chi-sq based on only polymorphic CpGs - BEST TEST since it controls for higher CpG mutation rate vs. rest of genome
print("####################################### chi-sq based on outlier SNPs in CpGs / nSNPs in CpGs #######################################")
#nDMLSNP / nCpGSNP * nCpGoutlier / nCpGSNP ?

nDMLSNP <- as.numeric(numOverlaps(DML, allSNP))
print(paste("number of DMLs in all SNPs (not just outliers)", nDMLSNP))

nCpGSNP <- as.numeric(numOverlaps(CPGwSNPcov, allSNP))
print(paste("number of SNPs in CpG sites covered by bs data (not just outliers)", nCpGSNP))

outlierCPGSNPwcov <- as.numeric(numOverlaps(CPGwSNPcov, outlierSNP))
print(paste("number of SNPs in outlier CpG sites covered by bs data", outlierCPGSNPwcov))

exp_overlap_poly <- (nDMLSNP / nCpGSNP) * (outlierCPGSNPwcov / nCpGSNP)
print(paste("expected DML-SNP overlap rate (polymorphic CpGs only)", exp_overlap_poly))

obs_overlap_poly <- overlaps / nCpGSNP
print(paste("observed DML-SNP overlap rate (polymorphic CpGs only)", obs_overlap_poly))

#nCpGSNP2 <- as.numeric(nCpGSNP * nCpGSNP)
nDMLSNP_nCpGoutlier <- as.numeric(nDMLSNP*outlierCPGSNPwcov)

#cont_tab<-rbind(c(nCpGSNP2, nDMLSNP_nCpGoutlier),
#                  c(nCpGSNP,overlaps))
#print(paste("contingency table", cont_tab))

#x2 <- chisq.test(cont_tab, correct=FALSE)
#print(paste("chi-sq results:"))
#print(x2)




nDMLSNP_nCpGoutlier <- as.numeric((nDMLSNP*outlierCPGSNPwcov)/nCpGSNP)

cont_tab<-rbind(c(nCpGSNP, nDMLSNP_nCpGoutlier),
                  c(nCpGSNP,overlaps))
print(paste("contingency table", cont_tab))

x2 <- chisq.test(cont_tab, correct=FALSE)
print(paste("chi-sq results:"))
print(x2)
