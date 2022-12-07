#!/bin/bash
#12_get_avg_DML_meth.sh
# get percent methylation per individual and per group from methylation matrix

METH_INPUT="/project/lbernatchez/users/clven7/whitefish/transcriptome/cliff/04_meth"
OUTPUT="10_DML_SNP_overlap_allelefreqs/allCpGs_"
METH_D="cliff_pct_meth_matrix_pct_meth_dwarf_means.txt"
OVERLAPS_D="10_DML_SNP_overlap_allelefreqs/allCpGs_CD_maf0.05_pctind0.75_maxdepth25_freqCsGs.bed"
METH_N="cliff_pct_meth_matrix_pct_meth_normal_means.txt"
OVERLAPS_N="10_DML_SNP_overlap_allelefreqs/allCpGs_CN_maf0.05_pctind0.75_maxdepth25_freqCsGs.bed"


module load bedtools

base_d="$(ls -1 "$METH_INPUT"/"$METH_D" | cut -d "_" -f2,3,4)"
awk 'NR>1 {print}' "$METH_INPUT"/"$METH_D" |
bedtools intersect -a stdin -b "$OVERLAPS_D" -nonamecheck -wo |
awk 'NR == 1 {print "chrom", "start", "end", "pct_meth_dwarf", "freq_C_dwarf"} 
	NR > 1 {print $1, $2, $3, $4, $11}' > "$OUTPUT"$(basename $base_d)"_pct_meth_freqC_DMLs_dwarf.txt"

base_n="$(ls -1 "$METH_INPUT"/"$METH_N" | cut -d "_" -f2,3,4)"
awk 'NR>1 {print}' "$METH_INPUT"/"$METH_N" |
bedtools intersect -a stdin -b "$OVERLAPS_N" -nonamecheck -wo |
awk 'NR == 1 {print "chrom", "start", "end", "pct_meth_normal", "freq_C_normal"} 
	NR > 1 {print $1, $2, $3, $4, $11}' > "$OUTPUT"$(basename $base_n)"_pct_meth_freqC_DMLs_normal.txt"