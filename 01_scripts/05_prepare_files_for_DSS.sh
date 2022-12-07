#!/bin/bash
# 02_prepare_files_for_DSS.sh
#srun -c 1 --mem 10G -p medium --time 7-00:00:00 -J prep_DSS -o 05_prep_DSS_%j.log ./01_scripts/05_prepare_files_for_DSS.sh &

INPUT="06_fully_filtered_bedGraphs"
OUTPUT="07_DSS_input"

for i in $(ls -1 "$INPUT"/*.bedGraph.gz | cut -d "_" -f4) 
do
    echo $(basename $i)
    gunzip -c "$INPUT"/"$(basename $i)"_fully_filtered.bedGraph.gz | 
#   awk -v FS="\t" -v OFS="\t" 'NR {print $1, $2, $7, $5}' |
	  awk -v FS="\t" -v OFS="\t" 'BEGIN{print "chr", "pos", "N", "X"}; NR {print $1, $2, $7, $5}' |
    gzip -c - > "$OUTPUT"/"$(basename $i)"".txt.gz"
done