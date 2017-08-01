#!/bin/bash
## Author S. Monzon
## version v2.0

# Test whether the script is being executed with sge or not.
if [ -z $SGE_TASK_ID ]; then
   	use_sge=0
else
   	use_sge=1
fi


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
set -x

## Usage

if [ $# != 13 -a "$use_sge" == "1" ]; then
   	echo "usage: ............"
   	exit
fi

#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed
set -x
echo `date`

# Variables

DIR_BAM=$1
THREADS=$2
REF_PATH=$3
OUTPUT_DIR=$4
SAMPLE_NAMES=$5
BAM_NAMES=$6
OUTPUT_VCF_NAME=$7
OUTPUT_SNPS_NAME=$8
OUTPUT_SNPS_NAME_FIL=$9
OUTPUT_INDELS_NAME=${10}
OUTPUT_INDELS_NAME_FIL=${11}
OUTPUT_VCF_FIL_NAME=${12}
GATK_PATH=${13}


mkdir -p $OUTPUT_DIR/variants
echo $BAM_NAMES | tr ":" "\n" | awk -v prefix=$DIR_BAM '{print prefix "/" $0}' > $OUTPUT_DIR/bam.list

java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
      -T HaplotypeCaller \
      -R $REF_PATH \
      -I $OUTPUT_DIR/bam.list \
      -o $OUTPUT_DIR/variants/$OUTPUT_VCF_NAME \
      -stand_call_conf 30.0 \
      -stand_emit_conf 10.0 \
      -ploidy 1 \
      -S LENIENT \
      -log $OUTPUT_DIR/$OUTPUT_VCF_NAME-HaplotypeCaller.log

java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
      -R $REF_PATH \
      -T SelectVariants \
      -V $OUTPUT_DIR/variants/$OUTPUT_VCF_NAME \
      -o $OUTPUT_DIR/variants/$OUTPUT_SNPS_NAME \
      -selectType SNP \
      -nt $THREADS \
      -S LENIENT \
      -log $OUTPUT_DIR/$OUTPUT_VCF_NAME-selectSNP.log


java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
	-R $REF_PATH \
	-T VariantFiltration \
	-V $OUTPUT_DIR/variants/$OUTPUT_SNPS_NAME \
	-o $OUTPUT_DIR/variants/$OUTPUT_SNPS_NAME_FIL \
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
	-log $OUTPUT_DIR/$OUTPUT_VCF_NAME-filterSNPs.log

echo -e "Select and Filter Indels"
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
	-R $REF_PATH \
	-T SelectVariants \
	-V $OUTPUT_DIR/variants/$OUTPUT_VCF_NAME \
	-o $OUTPUT_DIR/variants/$OUTPUT_INDELS_NAME \
	-selectType INDEL \
	-nt $THREADS \
	-S LENIENT \
	-log $OUTPUT_DIR/$OUTPUT_VCF_NAME-selectIndels.log

java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R $REF_PATH \
	-V $OUTPUT_DIR/variants/$OUTPUT_INDELS_NAME \
	-o $OUTPUT_DIR/variants/$OUTPUT_INDELS_NAME_FIL \
	--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0" \
    --filterName "IndelFilters" \
	-S LENIENT \
	-log $OUTPUT_DIR/$OUTPUT_VCF_NAME-filterIndels.log

echo -e "Combine snps and indels vcf"
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
    -R $REF_PATH \
    -T  CombineVariants \
    --genotypemergeoption UNSORTED \
    --variant $OUTPUT_DIR/variants/$OUTPUT_SNPS_NAME_FIL \
    --variant $OUTPUT_DIR/variants/$OUTPUT_INDELS_NAME_FIL \
    -o $OUTPUT_DIR/variants/$OUTPUT_VCF_FIL_NAME \
    -log $OUTPUT_DIR/$OUTPUT_VCF_NAME-CombineVCF.log
