#!/usr/bin/env Rscript
#conda activate R403b
#srun -c 1 --mem 30G -p medium --time 7-00:00:00 -J DSS -o 06_DSS_%j.log Rscript ./01_scripts/06_DSS_model_no_filtration.R &

library(DSS)
library(data.table)
library(regioneR)
library(dplyr)

# load files
path = "07_DSS_input/"
files = list.files(path=path, pattern=paste0("*.txt.gz"))
make_df = lapply(paste0(path, files), function (x) fread(x, header=T))

# make BSseq object
#BSobj=makeBSseqData(dat=make_df, sampleNames=files)
BSobj = makeBSseqData(make_df, basename(files))
BSobj

save.image()

# experimental design
sample_info <- read.table("whitefish_cliff_with_vars_whitelist_DSS.txt", header=T)
vars <- sample_info[,1:3]

# DML test
DMLfit = DMLfit.multiFactor(BSobj, design=vars, formula=~ecotype)

colnames(DMLfit$X)
DMLtest.ecotype = DMLtest.multiFactor(DMLfit, coef="ecotypenormal")

# if you want to sort DMLs
#ix=sort(DMLtest.ecotype[,"pvals"], index.return=TRUE)$ix
#head(DMLtest.ecotype[ix,])

# if you want a list of all CpG sites analyzed and data on them (for chisq)
write.table(DMLtest.ecotype, file="08_DSS_results/cliff_DMLs_all_CpGs.txt", sep="\t", quote=FALSE)

# filter DMLs by p-value and write table
sig_DMLs <- DMLtest.ecotype[which(DMLtest.ecotype$fdrs < 0.05),]
write.table(sig_DMLs, file="08_DSS_results/cliff_DMLs_p0.05.txt", sep="\t", quote=FALSE)

DMRtest = callDMR(DMLtest.ecotype, p.threshold=0.001)
write.table(DMRtest, file="08_DSS_results/cliff_DMRs_p0.05.txt", sep="\t", quote=FALSE)

# get heatmap betas
dmrs <- toGRanges(DMRtest)

### to get DMR beta values:
### from Kyle's script: https://github.com/kylewellband/methylUtil/blob/master/01_scripts/04_DML_DMR_plotting.R
    hits <- findOverlaps(BSobj, dmrs, ignore.strand = TRUE)
    dmr_cov <- aggregate(getCoverage(BSobj, type = "Cov")[queryHits(hits),], by = list(subjectHits(hits)), FUN = sum, na.rm = TRUE)[,-1]
    dmr_M <- aggregate(getCoverage(BSobj, type = "M")[queryHits(hits),], by = list(subjectHits(hits)), FUN = sum, na.rm = TRUE)[,-1]
    mcols(dmrs) <- cbind(mcols(dmrs), dmr_M / dmr_cov)
#	fwrite(as.data.frame(dmrs), "08_DSS_results/cliff_DMRs_p0.05_heatmap_mean_Betas.txt", sep = "\t")
 
dat <- as.data.frame(dmrs)
datf <- dat %>% mutate_all(na_if,"")

write.table(datf, "08_DSS_results/cliff_DMRs_p0.05_heatmap_mean_Betas.txt", sep = "\t", quote=F, row.names=F)