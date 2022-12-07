#!/bin/bash
# bedGraph_coverage_filtration.sh
#srun -c 1 --mem 10G -p medium --time 7-00:00:00 -J cov_filt -o 01_cov_filt_%j.log ./01_scripts/02_bedGraph_coverage_filtration.sh &

INPUT="/project/lbernatchez/users/clven7/whitefish/am/ga_methyl_cliff/04_whitelisted_bedGraphs"
OUTPUT="/project/lbernatchez/users/clven7/whitefish/am/ga_methyl_cliff/05_coverage_filtered_bedGraphs"

for i in $(ls -1 "$INPUT"/*whitelisted.bedGraph.gz | cut -d "_" -f5) 
do
    echo $(basename $i)
    gunzip -c "$INPUT"/"$(basename $i)"_SNPsexcluded2_withscaffolds_whitelisted.bedGraph.gz | 
		awk -v FS="\t" -v OFS="\t" 'NR {print $1, $2, $3, $4, $5, $6, ($5+$6)}' |
		awk -v FS="\t" -v OFS="\t" '{ if (4<$7 && 101>$7) {print} } ' |
#		awk -v FS="\t" -v OFS="\t" '{ if (4<$7) {print} } ' |
#    gzip -c - > "$OUTPUT"/"$(basename $i)""_SNPsexcluded2_noscaffolds_coveragefiltered.bedGraph.gz"
    gzip -c - > "$OUTPUT"/"$(basename $i)""_SNPsexcluded2_withscaffolds_coveragefiltered_min5_max100.bedGraph.gz"
done




#awk -v FS="\t" -v OFS="\t" 'NR {print $1, $2, $3, $4, $5, $6, ($5+$6)}' CD17_SNPsexcluded2_noscaffolds_whitelisted.bedGraph > CD17_SNPsexcluded2_noscaffolds_whitelisted_testsum.bedGraph
#awk -v FS="\t" -v OFS="\t" '{ if (4<$7 && $7<101) {print} } '  CD17_SNPsexcluded2_noscaffolds_whitelisted_testsum.bedGraph > CD17_SNPsexcluded2_noscaffolds_whitelisted_testsum_covfilter.bedGraph
