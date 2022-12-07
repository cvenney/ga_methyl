#!/bin/bash
# 01_whitelist_CpGs_with_SNPs.sh
#srun -c 1 --mem 15G -p medium --time 7-00:00:00 -J whitelist -o whitelist_%j.log ./01_scripts/01_whitelist_CpGs_with_SNPs.sh &

# no scaffolds
#INPUT="/project/lbernatchez/users/clven7/whitefish/am/methylutil-master_whitelist/03_raw_bedGraphs"

# with scaffolds
INPUT="/project/lbernatchez/users/clven7/whitefish/am/bwa-meth_pipeline/07_methyl_dackel_SNPsincluded2"
WHITELIST="/project/lbernatchez/users/camer78/coregonus/Cliff_ANGSD/03_saf_maf_gl_all/CpG_filtering"
OUTPUT="/project/lbernatchez/users/clven7/whitefish/am/ga_methyl_cliff/04_whitelisted_bedGraphs"

module load bedtools

for i in $(ls -1 "$WHITELIST"/*_min0.7_FULLwhitelist.bed | cut -d "_" -f7) 
do
	echo $(basename $i)
#	gunzip -k "$INPUT"/"$(basename $i)"_SNPsexcluded2_noscaffolds.bedGraph.gz |
  bedtools intersect -nonamecheck -a "$INPUT"/HI*"$(basename $i)".dedup_CpG_merged.bedGraph.gz -b "$WHITELIST"/"$(basename $i)"_min0.7_FULLwhitelist.bed -f 1.0     -wa |
    gzip -c - > "$OUTPUT"/"$(basename $i)""_SNPsexcluded2_withscaffolds_whitelisted.bedGraph.gz"
done