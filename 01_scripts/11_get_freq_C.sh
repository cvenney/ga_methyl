#!/bin/bash
# 11_get_freq_C_v2.sh
# convert maf to bed format and restrict to outlier SNP-DML overlaps, then get freq C for each (0-1)
# 10 G memory is sufficient

INPUT_D="/project/lbernatchez/users/camer78/coregonus/Cliff_ANGSD/06_saf_maf_by_pop/CD"
MAF_D="CD_maf0.05_pctind0.75_maxdepth25.mafs"
INPUT_N="/project/lbernatchez/users/camer78/coregonus/Cliff_ANGSD/06_saf_maf_by_pop/CN"
MAF_N="CN_maf0.05_pctind0.75_maxdepth25.mafs"

Cs="/project/lbernatchez/users/clven7/whitefish/am/ga_permutation_cliff/02_genome/all_CpG_3_columns_Cs.bed"
Gs="/project/lbernatchez/users/clven7/whitefish/am/ga_permutation_cliff/02_genome/all_CpG_3_columns_Gs.bed"

OVERLAPS="09_DML_SNP_overlap/overlaps_cliff_DMLs_p0.05_SNPs.txt"
OUTPUT="10_DML_SNP_overlap_allelefreqs"

module load bedtools

# need separate lists of all Cs and all Gs in CpG sites in genome (for relative position)
# North America
#awk -v FS="\t" -v OFS="\t" 'NR {print $1, $2, $3-1}' "/project/lbernatchez/users/clven7/whitefish/am/ga_permutation_cliff/02_genome/all_CpG_3_columns.bed" > "/project/lbernatchez/users/clven7/whitefish/am/ga_permutation_cliff/02_genome/all_CpG_3_columns_Cs.bed"
#awk -v FS="\t" -v OFS="\t" 'NR {print $1, $2+1, $3}' "/project/lbernatchez/users/clven7/whitefish/am/ga_permutation_cliff/02_genome/all_CpG_3_columns.bed" > "/project/lbernatchez/users/clven7/whitefish/am/ga_permutation_cliff/02_genome/all_CpG_3_columns_Gs.bed"

### dwarf
# Cs
base="$(ls -1 "$INPUT_D"/"$MAF_D" | cut -d "." -f1,2,3)"
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2-1, $2, $3, $4, $5, $6, $7}' $base.mafs |
bedtools intersect -a stdin -b "$Cs" |
#bedtools intersect -a stdin -b "$OVERLAPS" > "$OUTPUT"/$(basename $base)"_overlaps.bed"
bedtools intersect -a stdin -b "$OVERLAPS" |
awk -v FS="\t" -v OFS="\t" '{if ($4 == "C") print $1, $2, $3, $4, $5, $6, 1-$7;
	else print $1, $2, $3, $4, $5, $6, $7;}' - > "$OUTPUT"/$(basename $base)"_freqC.bed"
	
# Gs
base="$(ls -1 "$INPUT_D"/"$MAF_D" | cut -d "." -f1,2,3)"
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2-1, $2, $3, $4, $5, $6, $7}' $base.mafs |
bedtools intersect -a stdin -b "$Gs" |
#bedtools intersect -a stdin -b "$OVERLAPS" > "$OUTPUT"/$(basename $base)"_overlaps.bed"
bedtools intersect -a stdin -b "$OVERLAPS" |
awk -v FS="\t" -v OFS="\t" '{if ($4 == "G") print $1, $2, $3, $4, $5, $6, 1-$7;
	else print $1, $2, $3, $4, $5, $6, $7;}' - > "$OUTPUT"/$(basename $base)"_freqG.bed"


### normal
# Cs
base="$(ls -1 "$INPUT_N"/"$MAF_N" | cut -d "." -f1,2,3)"
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2-1, $2, $3, $4, $5, $6, $7}' $base.mafs |
bedtools intersect -a stdin -b "$Cs" |
#bedtools intersect -a stdin -b "$OVERLAPS" > "$OUTPUT"/$(basename $base)"_overlaps.bed"
bedtools intersect -a stdin -b "$OVERLAPS" |
awk -v FS="\t" -v OFS="\t" '{if ($4 == "C") print $1, $2, $3, $4, $5, $6, 1-$7;
	else print $1, $2, $3, $4, $5, $6, $7;}' - > "$OUTPUT"/$(basename $base)"_freqC.bed"
	
# Gs
base="$(ls -1 "$INPUT_N"/"$MAF_N" | cut -d "." -f1,2,3)"
awk -v FS="\t" -v OFS="\t" 'NR>1 {print $1, $2-1, $2, $3, $4, $5, $6, $7}' $base.mafs |
bedtools intersect -a stdin -b "$Gs" |
#bedtools intersect -a stdin -b "$OVERLAPS" > "$OUTPUT"/$(basename $base)"_overlaps.bed"
bedtools intersect -a stdin -b "$OVERLAPS" |
awk -v FS="\t" -v OFS="\t" '{if ($4 == "G") print $1, $2, $3, $4, $5, $6, 1-$7;
	else print $1, $2, $3, $4, $5, $6, $7;}' - > "$OUTPUT"/$(basename $base)"_freqG.bed"
	
	
cat "$OUTPUT"/CD_*freqC.bed "$OUTPUT"/CD_*freqG.bed > "$OUTPUT"/CD_maf0.05_pctind0.75_maxdepth25_freqCsGs.bed
cat "$OUTPUT"/CN_*freqC.bed "$OUTPUT"/CN_*freqG.bed > "$OUTPUT"/CN_maf0.05_pctind0.75_maxdepth25_freqCsGs.bed