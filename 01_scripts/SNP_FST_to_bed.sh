#SNP_FST_to_bed.sh

INPUT="02_SNPs/CN_CD_maf0.05_pctind0.75_maxdepth25_bypos_5pct.fst"
OUTPUT="02_SNPs/CN_CD_maf0.05_pctind0.75_maxdepth25_bypos_5pct.bed"

awk 'NR==1 {print "#chrom", "start", "stop", "fst"} NR > 1 {print $2, $3-1, $3, $6}' "$INPUT" > "$OUTPUT"