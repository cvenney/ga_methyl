#!/bin/bash
#07_whitefish_overlap_DMLs_and_SNPs.sh
#srun -c 1 --mem 10G -p small --time 1-00:00:00 -J 07_overlap -o 07_overlap_%j.log ./01_scripts/07_whitefish_overlap_DMLs_and_SNPs.sh &

DSS_RESULTS="08_DSS_results"
SNPs="02_SNPs/CN_CD_maf0.05_pctind0.75_maxdepth25_bypos_5pct.bed"

module load bedtools 

### DMLs
for i in $(ls -1 "$DSS_RESULTS"/*_DMLs_*.txt | sed 's/.txt//g')


do
    echo $(basename $i)

	# convert DML file to bed format - DSS files are based on start-1 from bed format
	awk -v FS="\t" -v OFS="\t" 'NR>1 {print $2, $3, $3+2}' "$DSS_RESULTS"/$(basename $i).txt > "$DSS_RESULTS"/$(basename $i).bed

	# intersect file with SNPs
	bedtools intersect -a "$DSS_RESULTS"/$(basename $i).bed -b "$SNPs" -nonamecheck > 09_DML_SNP_overlap/overlaps_"$(basename $i)"_"SNPs.txt"

done

### DMRs
for i in $(ls -1 "$DSS_RESULTS"/*_DMRs_*.txt | sed 's/.txt//g')


do
    echo $(basename $i)

	# convert DMR file to bed format
	awk -v FS="\t" -v OFS="\t" 'NR>1 {print $2, $3, $4+2}' "$DSS_RESULTS"/$(basename $i).txt > "$DSS_RESULTS"/$(basename $i).bed

	# intersect file with SNPs
	bedtools intersect -a "$DSS_RESULTS"/$(basename $i).bed -b "$SNPs" -nonamecheck > 09_DML_SNP_overlap/overlaps_"$(basename $i)"_"SNPs.txt"

done


### full list of CpG sites (for DMR-SNP chisq test)
awk 'NR > 1 {print $2, $3, $3+2}' OFS='\t' 08_DSS_results/cliff_DMLs_all_CpGs.txt > 08_DSS_results/cliff_DMLs_all_CpGs.bed