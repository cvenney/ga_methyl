#!/bin/bash
#20_overlap_genes_with_sig_markers.sh

GENES="/project/lbernatchez/users/clven7/whitefish/am/genome_annotation_table_simplified_na_noheader.txt"
DML="08_DSS_results/cliff_DMLs_p0.05.bed"
OUTLIERS="02_SNPs/CN_CD_maf0.05_pctind0.75_maxdepth25_bypos_5pct.bed"
DML_SNP="09_DML_SNP_overlap/overlaps_cliff_DMLs_p0.05_SNPs.txt"
DMR="08_DSS_results/cliff_DMRs_p0.05.bed"
OUTPUT="12_gene_associations"

# need to first reformat gene table for use in bedtools
#awk -v FS="\t" -v OFS="\t" 'NR > 1 {print}' /project/lbernatchez/users/clven7/whitefish/am/genome_annotation_table_simplified_na.tsv > "$GENES"

module load bedtools

bedtools intersect -a "$GENES" -b "$DML" -wb -wa |
awk -v FS="\t" -v OFS="\t" '{print $8}' | sort | uniq |
sed 's/ /_/g' | sed "s/[^[:alnum:]_-]//g" > "$OUTPUT"/cliff_DML_genes_unique.txt

bedtools intersect -a "$GENES" -b "$DMR" -wb -wa |
awk -v FS="\t" -v OFS="\t" '{print $8}' | sort | uniq |
sed 's/ /_/g' | sed "s/[^[:alnum:]_-]//g" > "$OUTPUT"/cliff_DMR_genes_unique.txt

bedtools intersect -a "$GENES" -b "$DML_SNP" -wb -wa |
awk -v FS="\t" -v OFS="\t" '{print $8}' | sort | uniq |
sed 's/ /_/g' | sed "s/[^[:alnum:]_-]//g" > "$OUTPUT"/cliff_DML-SNP_overlaps_genes_unique.txt

bedtools intersect -a "$GENES" -b "$OUTLIERS" -wb -wa |
awk -v FS="\t" -v OFS="\t" '{print $8}' | sort | uniq |
sed 's/ /_/g' | sed "s/[^[:alnum:]_-]//g" > "$OUTPUT"/cliff_outlier_SNP_overlaps_genes_unique.txt