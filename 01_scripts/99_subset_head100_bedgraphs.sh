#!/bin/bash
# subset_head100_bedgraphs.sh
INPUT="/project/lbernatchez/users/clven7/whitefish/am/methylutil-master_whitelist/05_coverage_filtered_bedGraphs"

for i in $(ls -1 "$INPUT"/*coveragefiltered.bedGraph.gz | cut -d "_" -f5) 
do
#    echo $(basename $i)
	gunzip -c "$INPUT"/"$(basename $i)"_SNPsexcluded2_noscaffolds_coveragefiltered.bedGraph.gz  | 
	head -2000 - > "$INPUT"/"$(basename $i)"_SNPsexcluded2_noscaffolds_coveragefiltered_head200.bedGraph
done