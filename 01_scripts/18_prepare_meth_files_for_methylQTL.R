#!/usr/bin/env Rscript
#18_prepare_meth_files_for_methylQTL.R
# conda activate R42
# srun -c 16 --mem 40G --time 7-00:00:00 -p medium -o 18_GEM_methylQTL_%j.log Rscript ./01_scripts/18_GEM_methylQTL.R &

# first format methylation matrix into frequency

full <- read.table("11_methylQTL/cliff_pct_meth_matrix_full.txt", header=T)
dml <- read.table("11_methylQTL/cliff_DMLs_p0.05_positions.txt", header=T)

library(dplyr)

meth2 <- full
meth2[2] = meth2[2]/100
meth2[3] = meth2[3]/100
meth2[4] = meth2[4]/100
meth2[5] = meth2[5]/100
meth2[6] = meth2[6]/100
meth2[7] = meth2[7]/100
meth2[8] = meth2[8]/100
meth2[9] = meth2[9]/100
meth2[10] = meth2[10]/100
meth2[11] = meth2[11]/100
meth2[12] = meth2[12]/100
meth2[13] = meth2[13]/100
meth2[14] = meth2[14]/100
meth2[15] = meth2[15]/100
meth2[16] = meth2[16]/100
meth2[17] = meth2[17]/100

#write.table(meth2, "11_methylQTL/cliff_pct_meth_matrix_full_freq.txt", sep = "\t", row.names=F)

# extract info on DMLs only

dml_freq <- inner_join(meth2, dml)

write.table(dml_freq, "11_methylQTL/cliff_pct_meth_matrix_dml_freq.txt", sep = "\t", row.names=F, quote=F)