#!/bin/bash
set -e
set -u
set -o  pipefail

usage(){
cat << EOF
Usage: 
	bash $(echo $0) [-r <RSCRIPT>][-d <DUPLICATION_FILE>] [-c <HICHIP_FILE>][-t <TSS_FILE>][-g <GENEBODY_FILE>]
       -r Rscript	Rscript file [required].
       -d DUPLICATION_FILE	Duplication or region file [required].
       -c HICHIP_FILE	HiChIP file from Hichipper pipeline or other bedpe file [required].
       -t TSS_FILE	TSS file for protein coding genes or other TSS file [required].
       -g GENEBODY_FILE	Gene body file for protein coding genes [required].
       -s GENE_ENH	Gene or region in Duplication to search; ALL gene and region by default. e.g.,"C4BPB","chr1:201970000-202085000" or "C4BPB:chr1:201970000-202085000" [optional]

Example:
	bash ./src/RENC_v2.sh -r ./src/RENC_v2.R  
	-d ./example/Duplication_STAD_chr1.bed  
	-c ./example/AGS_STAD_chr1.bedpe  
	-t ./reference/hg19/RefSeq_proteinCoding.tss.bed  
	-g ./reference/hg19/RefSeq_proteinCoding.body.bed 

EOF
}
GENE_ENH=""
while getopts ":hr:d:c:t:g:s:" arg
do
	case "$arg" in
        "r")
		    RSCRIPT="$OPTARG"
			;;
		"d")
		    DUPLICATION_FILE="$OPTARG"
			;;
		"c")
		    HICHIP_FILE="$OPTARG"
		  ;;
		"t")
			TSS_FILE="$OPTARG"
			;;
        "g")
			GENEBODY_FILE="$OPTARG"
			;;
        "s")
			GENE_ENH="$OPTARG"
			;;   
		"?")
			echo "Unkown option: $OPTARG"
			exit 1
			;;
		":")
			echo "No argument value for option $OPTARG"
			;;
		h)
			usage
			exit 0
			;;
		*)
			echo "Unknown error while processing options"
			usage
			exit 1
			;;
	esac
done

if [ $# != 12 ] && [ $# != 10 ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

if [ -z "${GENE_ENH}" ]; then
	GENE_ENH="ALL_RESULT"
fi

master_time_start=`date +%s` 
sample=${HICHIP_FILE##*/}
sample=${sample%%.*}
printf "\n---------Begin------------\n"
#echo ${sample}
bedtools pairtobed -a ${HICHIP_FILE} -b ${TSS_FILE} -type xor > ${sample}_tss.bed
## Non-gene duplication hotspots
bedtools intersect -a ${DUPLICATION_FILE} -b ${GENEBODY_FILE} -u -F 1.0 > ${sample}_gene.bed
bedtools subtract -a ${DUPLICATION_FILE} -b  ${sample}_gene.bed > ${sample}_non_gene.bed
bedtools pairtobed -a ${HICHIP_FILE} -b ${sample}_non_gene.bed > ${sample}_duplication.bedpe
bedtools pairtobed -a ${sample}_duplication.bedpe -b ${TSS_FILE} -type xor > ${sample}_non_gene_tss.bedpe
cat ${sample}_non_gene_tss.bedpe | awk '{if(($14>$5 && $15<$6 && $2>$10 && $2<$11)||($14>$5 && $15<$6 && $3>$10 && $3<$11)){print $0}}'  > ${sample}_non_gene_tss1.bedpe
cat ${sample}_non_gene_tss.bedpe | awk '{if(($14>$2 && $15<$3 && $5>$10 && $5<$11)||($14>$2 && $15<$3 && $6>$10 && $6<$11)){print $0}}'  > ${sample}_non_gene_tss2.bedpe
cat ${sample}_non_gene_tss1.bedpe ${sample}_non_gene_tss2.bedpe > ${sample}_no_gene.bedpe
rm -rf ${sample}_gene.bed ${sample}_non_gene.bed ${sample}_duplication.bedpe 
rm -rf ${sample}_non_gene_tss.bedpe ${sample}_non_gene_tss1.bedpe ${sample}_non_gene_tss2.bedpe

## Gene-harboring duplication hotspots
# TSS is outside duplication
bedtools intersect -a ${DUPLICATION_FILE} -b ${GENEBODY_FILE} -wa -wb -F 1.0 > ${sample}_gene.bed
bedtools pairtobed -a ${HICHIP_FILE} -b ${sample}_gene.bed > ${sample}_duplication.bedpe
bedtools pairtobed -a ${sample}_duplication.bedpe -b ${TSS_FILE} -type xor > ${sample}_gene_tss.bedpe
cat ${sample}_gene_tss.bedpe | awk '{if(($19!=$24 && $21>$5 && $22<$6 && $2>$10 && $2<$11)||($19!=$24 && $21>$5 && $22<$6 && $3>$10 && $3<$11)){print $0}}' > ${sample}_gene_tss1.bedpe
cat ${sample}_gene_tss.bedpe | awk '{if(($19!=$24 && $21>$2 && $22<$3 && $5>$10 && $5<$11)||($19!=$24 && $21>$2 && $22<$3 && $6>$10 && $6<$11)){print $0}}' > ${sample}_gene_tss2.bedpe
cat ${sample}_gene_tss1.bedpe ${sample}_gene_tss2.bedpe > ${sample}_gene_tss_out_dup.bedpe
# TSS is inside duplication
cat ${sample}_gene_tss.bedpe | awk '{if($19==$24 && $21>$5 && $22<$6 && $3>$10 && $3<$11){print $0}}' > ${sample}_gene_tss1.bedpe
cat ${sample}_gene_tss.bedpe | awk '{if($19==$24 && $21>$2 && $22<$3 && $5>$10 && $5<$11){print $0}}' > ${sample}_gene_tss2.bedpe
cat ${sample}_gene_tss1.bedpe ${sample}_gene_tss2.bedpe >  ${sample}_gene_tss_in_dup.bedpe
rm -rf ${sample}_gene.bed ${sample}_duplication.bedpe ${sample}_gene_tss.bedpe
rm -rf ${sample}_gene_tss1.bedpe ${sample}_gene_tss2.bedpe

mkdir -p output/bedpe_tss output/no_gene output/gene_tss_out_dup output/gene_tss_in_dup
mv *_tss.bed output/bedpe_tss
mv *_no_gene.bedpe output/no_gene
mv *_gene_tss_out_dup.bedpe output/gene_tss_out_dup
mv *_gene_tss_in_dup.bedpe output/gene_tss_in_dup

# RENC result
#echo "=====Rscript begin====="
echo "Run script: ${RSCRIPT} ${TSS_FILE} ${GENE_ENH}"
Rscript ${RSCRIPT} ${TSS_FILE} ${GENE_ENH}
rm -rf output

master_time_end=`date +%s`
master_time_exec=`echo "$master_time_end - $master_time_start"|bc`
master_time_exec=`echo "$master_time_exec"|bc`
echo Elapsed $master_time_exec sec
printf "\n---------Done-------------\n"

# Print search result
if [ ${GENE_ENH} != "ALL_RESULT" ]; then
	cat *search.txt
fi






