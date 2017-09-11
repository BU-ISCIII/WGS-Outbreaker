#!/bin/bash
## Author S. Monzon
## version v2.0

# Help
# usage: run_variantCalling.sh ....

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#Execute processing_config.sh
source $SCRIPTS_DIR/processing_config.sh

# VARIABLES

#USE_SGE=$1
#VARIANT_CALLING=$2
#DUPLICATES=$3
#OUTPUT_DIR=$4
#REF_PATH=$5
#THREADS=$6
#SAMPLE_NAMES=$7
#OUTPUT_BAM_NAMES=$8
#VCF_NAMES=$9
#sample_number=${10}
#OUTPUT_BAM_REALIGNED_NAMES=${11}
#OUTPUT_BAM_RECALIBRATED_NAMES=${12}
#GATK_PATH=${13}
#SNPS_NAME=${14}
#SNPS_FIL_NAME=${15}
#INDELS_NAME=${16}
#INDELS_FIL_NAME=${17}
#VCF_FIL_NAME=${18}
#KNOWN_SNPS=${19}
#KNOWN_INDELS=${20}

## Folder creation
echo -e "Creating $output_dir/variant_calling"
mkdir -p $output_dir/variant_calling

if [ "$use_sge" = "1" -a $duplicate_filter == "YES" ]; then
 	jobid=$(cat $output_dir/logs/jobids.txt | grep -w "PICARD" | cut -d ':' -f2 )
 	precalling_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid"
else
 	jobid=$(cat $output_dir/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )
 	precalling_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid"
fi



if [ $duplicate_filter == "YES" ]; then
	input_list=$duplicateBamArray_list
else
	input_list=$mappingArray_sorted_list
fi

mkdir -p $output_dir/variant_calling/variants_gatk
precalling_cmd="$SCRIPTS_DIR/gatk_preprocessing.sh \
	$threads \
	$output_dir/Alignment/BAM \
	$output_dir/variant_calling/variants_gatk \
	$samples \
  	$know_snps \
	$know_indels \
	$input_list \
	$realignedBamArray_list \
	$recalibratedBamArray_list \
	$ref_path \
	$gatk_path"


if [ $know_snps == "NO" ];then
	calling_cmd="$SCRIPTS_DIR/gatk_haploid.sh \
	$threads \
	$output_dir/variant_calling/variants_gatk/realignment \
	$output_dir/variant_calling/variants_gatk \
	$samples \
	$realignedBamArray_list \
	$vcffilArray_list \
	$vcfArray_list \
	$vcfsnpsArray_list \
	$vcfsnpsfilArray_list \
	$vcfindelsArray_list \
	$vcfindelsfilArray_list \
	$ref_path \
	$gatk_path"

else
  	calling_cmd="$SCRIPTS_DIR/gatk_haploid.sh \
	$threads \
	$output_dir/variant_calling/variants_gatk/recalibration \
	$output_dir/variant_calling/variants_gatk \
	$samples \
	$recalibratedBamArray_list \
	$vcfsnpsArray_list \
        $vcfsnpsfilArray_list \
        $vcfindelsArray_list \
        $vcfindelsfilArray_list \
        $vcffilArray_list \
        $ref_path \
        $gatk_path"

fi

if [ $variant_calling == "YES" ]; then
		if [ "$use_sge" = "1" ]; then
             	precalling=$( qsub $precalling_args -t 1-$sample_count -N $JOBNAME.CALLING $precalling_cmd)
    		jobid_precalling=$( echo $precalling | cut -d ' ' -f3 | cut -d '.' -f1 )
    		calling_args="${SGE_ARGS} -hold_jid $jobid_precalling"
    		callling=$( qsub $calling_Args -N $JOBNAME.CALLING $calling_cmd)
       		jobid_calling=$( echo $calling | cut -d ' ' -f3 | cut -d '.' -f1 )
       		echo -e "Variant Calling:$jobid_precalling - $jobid_calling \n" >> $output_dir/logs/jobids.txt
		else
        	for count in `seq 1 $sample_count`
        	do
        		echo "Running variant calling on sample $count"
        		precalling=$($precalling_cmd $count)
       		done
        	calling=$($calling_cmd)
      fi
fi
