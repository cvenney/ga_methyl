#!/bin/bash
# 17_convert_geno_meth_for_GEM.sh

GENO="/project/lbernatchez/users/camer78/coregonus/Cliff_ANGSD/03_saf_maf_gl_all/all_maf0.05_pctind0.75_maxdepth25_cutOff0.7.geno.gz"
METH="/project/lbernatchez/users/clven7/whitefish/transcriptome/cliff/04_meth/cliff_pct_meth_matrix_full.txt"
DML="/project/lbernatchez/users/clven7/whitefish/am/ga_methyl_cliff/08_DSS_results/cliff_DMLs_p0.05.txt"
OUTPUT="11_methylQTL"

gunzip -c "$GENO" |
	awk -v FS="\t" -v OFS="\t" 'BEGIN{print "id", "CD17", "CD18", "CD19", "CD20", "CD21", "CD22", "CD28", "CD32", "CN10", "CN11", "CN12", "CN14", "CN15", "CN5", "CN6", "ID9"}; {print $1 "_" $2, $3+1, $4+1, $5+1, $6+1, $7+1, $8+1, $9+1, $10+1, $11+1, $12+1, $13+1, $14+1, $15+1, $16+1, $17+1, $18+1}' |
	awk -v FS="\t" -v OFS="\t" 'NR > 1 {for(x=2;x<=NF;x++){gsub("0","NA",$x)}}1' > "$OUTPUT"/all_maf0.05_pctind0.75_maxdepth25_cutOff0.7.geno.GEM
	
awk -v FS="\t" -v OFS="\t" 'NR == 1 {print "id", $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19}; NR > 1 {print $1 "_" $2+1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19}' "$METH" > "$OUTPUT"/cliff_pct_meth_matrix_full.txt

awk -v FS="\t" -v OFS="\t" 'NR == 1 {print "id"}; NR > 1 {print $2 "_" $3+1}' "$DML" > "$OUTPUT"/cliff_DMLs_p0.05_positions.txt