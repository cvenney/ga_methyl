#!/usr/bin/env Rscript
#conda activate R403b
#usage: Rscript ./01_scripts/make_CpG_blacklist_cov_across_individuals.R 4 sample_info/cliff_dwarf.txt sample_info/cliff_normal.txt 05_coverage_filtered_bedGraphs/03_make_CpG_blacklist_cov_across_individuals.R
#srun -c 1 --mem 40G -p medium --time 7-00:00:00 -J cl_make_blacklist -o 03_make_blacklist_%j.log Rscript ./01_scripts/03_make_CpG_blacklist_cov_across_individuals.R 4 sample_info/cliff_dwarf.txt sample_info/cliff_normal.txt 05_coverage_filtered_bedGraphs/cliff_blacklist.bed &

library(data.table)
library(dplyr)

argv <- commandArgs(T)
max_missing <- as.numeric(argv[1])
sample_file <- argv[2]
sample_file_n <- argv[3]
output <- argv[4]

#max_missing <- 4
#sample_file <- "sample_info/indian_dwarf.txt"
#sample_file_n <- "sample_info/indian_normal.txt"
#output <- "05_coverage_filtered_bedGraphs/indian_blacklist_v2test.bed"
#vector_ind<-c("CD17", "CD18", "CD19", "CD20", "CD21", "CD22", "CD28", "CD32")

ds <- read.table(sample_file, header=F)
ns <- read.table(sample_file_n, header=F)
all_samples <- rbind(ds, ns)

vector_ind<-all_samples[,1]
#read ind 1
bedgraph<-fread(paste0("05_coverage_filtered_bedGraphs/", vector_ind[1],"_SNPsexcluded2_withscaffolds_coveragefiltered_min5_max100.bedGraph.gz"))[,1:4]
colnames(bedgraph)<-c("chr","start","stop",vector_ind[1])

# loop from ind 2 to ind n to get full matrix
# need matrix of ALL individuals to account for samples with sufficient coverage in group 1 and **NO** coverage in group 2
# otherwise the site doesn't appear in group 2's blacklist and won't be masked in the dataset!
for (i in 2 : length(vector_ind))
{
print(vector_ind[i])
#bedgraph_i<-fread(paste0(vector_ind[i],"_SNPsexcluded2_noscaffolds_coveragefiltered.bedGraph.gz")
bedgraph_i<-fread(paste0("05_coverage_filtered_bedGraphs/", vector_ind[i],"_SNPsexcluded2_withscaffolds_coveragefiltered_min5_max100.bedGraph.gz"))
bedgraph_4<-bedgraph_i[,1:4]
colnames(bedgraph_4)<-c("chr","start","stop",vector_ind[i])
bedgraph<-full_join(bedgraph,bedgraph_4)
}

# split dwarf and normal
dw_df <- bedgraph %>% select(chr, start, stop, ds[,1])
no_df <- bedgraph %>% select(chr, start, stop, ns[,1])

### DWARF
count_na <- function(x) sum(is.na(x))

dw_df_filtered <- dw_df %>%
	mutate(count_na = apply(., 1, count_na)) %>%
	filter(count_na > max_missing)
dim(dw_df)
dim(dw_df_filtered)


### NORMAL
count_na_n <- function(x) sum(is.na(x))

no_df_filtered <- no_df %>%
	mutate(count_na_n = apply(., 1, count_na_n)) %>%
	filter(count_na_n > max_missing)
dim(no_df)
dim(no_df_filtered)


### COMBINE 

full <- rbind(dw_df_filtered[,1:3], no_df_filtered[,1:3])
full_final <- distinct(full)
dim(full)
dim(full_final)

write.table(full_final[,1:3], output, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
