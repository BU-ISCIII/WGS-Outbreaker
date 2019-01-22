#!/bin/bash
## Author S. Monzon
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to run GATK joint genotype\n"
        echo -e "Usage: gatk_jointvcf.sh\n \tthreads\n \toutput_dir\n \tsamples_list\n \tvcf_list\n \tvcf_fil_list\n \tvcf_snp_list\n \tvcf_snp_fils_list\n \tvcf_indel_list\n \tvcf_indel_fil_list\n \treference_path\n \tgatk_path\n \tgVCF_list\n \tvcf_PASS_list\n \tmax_snp\n \twindows_size\n"
        exit
fi


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
#Print commands and their arguments as they are executed.
set -x

#VARIABLE

threads=$1
output_dir=$2
samples=$3
vcfArray_list=$4
vcffilArray_list=$5
vcfsnpsArray_list=$6
vcfsnpsfilArray_list=$7
vcfindelsArray_list=$8
vcfindelsfilArray_list=$9
vcfsnpPassCluster=${10}
ref_path=${11}
gatk_path=${12}
haplotypeGVCF_list=${13}
vcfsnpPass=${14}
max_snp=${15}
window_size=${16}

echo $haplotypeGVCF_list | tr ":" "\n" | awk -v prefix=$output_dir/variants '{print prefix "/" $0}' > $output_dir/gvcf.list

# Merge gVCF
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
        -T GenotypeGVCFs \
        -R $ref_path \
        --variant $output_dir/gvcf.list \
        -o $output_dir/variants/$vcfArray_list \
	-log $output_dir/$vcfArray_list-jointGVCF.log

# Select snp variants
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
      -R $ref_path \
      -T SelectVariants \
      -V $output_dir/variants/$vcfArray_list \
      -o $output_dir/variants/$vcfsnpsArray_list \
      -selectType SNP \
      -S LENIENT \
      -log $output_dir/$vcfArray_list-selectSNP.log
#-nt $threads \

# Filter snp variants
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
        -R $ref_path \
        -T VariantFiltration \
        -V $output_dir/variants/$vcfsnpsArray_list \
        -o $output_dir/variants/$vcfsnpsfilArray_list \
	--clusterSize $max_snp\
        --clusterWindowSize $window_size \
        --filterExpression "MQ < 40.0" \
        --filterName "RMSMappingQuality" \
        --filterExpression "DP < 5" \
        --filterName "MaxDepth" \
        --filterExpression "QD < 2.0" \
        --filterName "LowQD" \
        --filterExpression "FS > 60.0" \
        --filterName "p-value StrandBias" \
        --filterExpression "SOR > 3.0" \
        --filterName "StandOddRatio" \
        -S LENIENT \
        -log $output_dir/$vcfArray_list-filterSNPs.log

#Select PASS and SnpCluster with vcftools

echo -e "Select PASS and SnpCluster"

vcftools --vcf $output_dir/variants/$vcfsnpsfilArray_list --remove-filtered 'RMSMappingQuality' --remove-filtered 'LowCoverage' --remove-filtered 'LowQD' --remove-filtered 'p-value StrandBias' --remove-filtered 'StandOddRatio' --recode -c > $output_dir/variants/$vcfsnpPassCluster

# Select only pass snp with vcftools
echo -e "Select PASS snp"

vcftools --vcf $output_dir/variants/$vcfsnpsfilArray_list --remove-filtered-all --recode -c > $output_dir/variants/$vcfsnpPass

# Select indels variants
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

# Filter indels vatiants
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
        -T VariantFiltration \
        -R $ref_path \
        -V $output_dir/variants/$vcfindelsArray_list \
        -o $output_dir/variants/$vcfindelsfilArray_list \
        --filterExpression "QD < 2.0 || FS > 200.0 || SOR > 10.0" \
    --filterName "IndelFilters" \
        -S LENIENT \
        -log $output_dir/$vcfArray_list-filterIndels.log

# Combine snp and indels in one vcf file
echo -e "Combine snps and indels vcf"
java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $gatk_path/GenomeAnalysisTK.jar \
    -R $ref_path \
    -T  CombineVariants \
    --genotypemergeoption UNSORTED \
    --variant $output_dir/variants/$vcfsnpsfilArray_list \
    --variant $output_dir/variants/$vcfindelsfilArray_list \
    -o $output_dir/variants/$vcffilArray_list \
    -log $output_dir/$vcfArray_list-CombineVCF.log

