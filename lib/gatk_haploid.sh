#!/bin/bash
## Author S. Monzon
## version v2.0


# Test whether the script is being executed with sge or not.
if [ -z $sge_task_id ]; then
   	use_sge=0
else
   	use_sge=1
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


## Usage

if [ $# != 13 -a "$use_sge" == "1" ]; then
   	echo "usage: ............"
   	exit
fi

echo `date`

# VARIABLES

threads=$1
input_dir=$2
output_dir=$3
samples=$4
BamArray_list=$5
vcfArray_list=$6
vcffilArray_list=$7
vcfsnpsArray_list=$8
vcfsnpsfilArray_list=$9
vcfindelsArray_list=${10}
vcfindelsfilArray_list=${11}
ref_path=${12}
gatk_path=${13}


mkdir -p $output_dir/variants
echo $BamArray_list | tr ":" "\n" | awk -v prefix=$input_dir '{print prefix "/" $0}' > $output_dir/bam.list

java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
      -T HaplotypeCaller \
      -R $ref_path \
      -I $output_dir/bam.list \
      -o $output_dir/variants/$vcfArray_list \
      -stand_call_conf 30.0 \
      -stand_emit_conf 10.0 \
      -ploidy 1 \
      -S LENIENT \
      -log $output_dir/$vcfArray_list-HaplotypeCaller.log

java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
      -R $ref_path \
      -T SelectVariants \
      -V $output_dir/variants/$vcfArray_list \
      -o $output_dir/variants/$vcfsnpsArray_list \
      -selectType SNP \
      -nt $threads \
      -S LENIENT \
      -log $output_dir/$vcfArray_list-selectSNP.log


java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
	-R $ref_path \
	-T VariantFiltration \
	-V $output_dir/variants/$vcfsnpsArray_list \
	-o $output_dir/variants/$vcfsnpsfilArray_list \
	--clusterWindowSize 10 \
	--filterExpression "MQ < 40" \
	--filterName "RMSMappingQuality" \
	--filterExpression "DP <5 " \
	--filterName "LowCoverage" \
	--filterExpression "QD <2.0 " \
	--filterName "LowQD" \
	--filterExpression "FS >60.0 " \
	--filterName "p-value StrandBias" \
	--filterExpression "MQRankSum < -12.5" \
	--filterName "MappingQualityRankSumTest" \
	--filterExpression "ReadPosRankSum < -8.0" \
	--filterName "VariantReadPosEnd" \
	--filterExpression "SOR > 4.0" \
	--filterName "StrandOddRank" \
	-S LENIENT \
	-log $output_dir/$vcfArray_list-filterSNPs.log

echo -e "Select and Filter Indels"
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
	-R $ref_path \
	-T SelectVariants \
	-V $output_dir/variants/$vcfArray_list \
	-o $output_dir/variants/$vcfindelsArray_list \
	-selectType INDEL \
	-nt $threads \
	-S LENIENT \
	-log $output_dir/$vcfArray_list-selectIndels.log

java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $ref_path \
	-V $output_dir/variants/$vcfindelsArray_list \
	-o $output_dir/variants/$vcfindelsfilArray_list \
	--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
    --filterName "IndelFilters" \
	-S LENIENT \
	-log $output_dir/$vcfArray_list-filterIndels.log

echo -e "Combine snps and indels vcf"
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
    -R $ref_path \
    -T  CombineVariants \
    --genotypemergeoption UNSORTED \
    --variant $output_dir/variants/$vcfsnpsfilArray_list \
    --variant $output_dir/variants/$vcfindelsfilArray_list \
    -o $output_dir/variants/$vcffilArray_list \
    -log $output_dir/$vcfArray_list-CombineVCF.log
