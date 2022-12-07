#!/bin/bash
# blacklist_CpGs_with_lov_num_ind.sh
#srun -c 1 --mem 20G -p medium --time 7-00:00:00 -J blacklist -o 04_blacklist_%j.log ./01_scripts/04_blacklist_CpGs_with_low_num_ind.sh &

BLACKLIST="05_coverage_filtered_bedGraphs/cliff_blacklist.bed"
INPUT="05_coverage_filtered_bedGraphs"
OUTPUT="06_fully_filtered_bedGraphs"

module load bedtools

for i in $(ls -1 "$INPUT"/*_SNPsexcluded2_withscaffolds_coveragefiltered_min5_max100.bedGraph.gz | cut -d "_" -f4) 

do
	echo $(basename $i)
	gunzip -c "$INPUT"/"$(basename $i)"_SNPsexcluded2_withscaffolds_coveragefiltered_min5_max100.bedGraph.gz  |
	bedtools intersect -nonamecheck -a stdin -b "$BLACKLIST" -wa -v |
    gzip -c - > "$OUTPUT"/"$(basename $i)""_fully_filtered.bedGraph.gz" 
done