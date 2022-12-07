# get beta values and plot heatmap
# assumes bsseq object is saved in R (from earlier scripts)
# if not, load in bsseq object as in 06_DSS

library(ComplexHeatmap)
library(circlize)
library(GenomicAlignments)
library(bsseq)
library(regioneR)
library(dplyr)

dmrs <- toGRanges(paste0("08_DSS_results/cliff_DMRs_p0.05.bed")

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