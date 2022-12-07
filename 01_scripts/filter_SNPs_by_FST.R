#filter_SNPs_by_FST.R

library(dplyr)

dat <- read.table("/project/lbernatchez/users/camer78/coregonus/Cliff_ANGSD/07_fst_by_pop_pair/CN_CD_maf0.05_pctind0.75_maxdepth25.bypos.fst", header=T)

sorted <- dat %>% arrange(desc(FST))

subset <- sorted[1:(nrow(sorted)*0.05),]

write.table(subset, "02_SNPs/CN_CD_maf0.05_pctind0.75_maxdepth25_bypos_5pct.fst")