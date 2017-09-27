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

vcf_to_tsv_cmd="$SCRIPTS_DIR/vcf_to_tsv.sh \
	$threads
	$output_dir/variant_calling/variants_gatk/variants \
	$vcfsnpsfilArray_list \
	$tsv_file"

tsv_to_msa_cmd="$SCRIPTS_DIR/tsv_to_msa.sh \
	$threads
        $output_dir/variant_calling/variants_gatk/variants \
	$tsv_file
	$msa_file"

if [ $variant_calling == "YES" ]; then
		if [ "$use_sge" = "1" ]; then
             	precalling=$( qsub $precalling_args -t 1-$sample_count -N $JOBNAME.CALLING $precalling_cmd)
    		jobid_precalling=$( echo $precalling | cut -d ' ' -f3 | cut -d '.' -f1 )
    		calling_args="${SGE_ARGS} -hold_jid $jobid_precalling"
    		calling=$( qsub $calling_Args -N $JOBNAME.CALLING $calling_cmd)
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

if [ $vcf_to_msa == "YES" ]; then
	if [ "$use_sge" = "1" ]; then
	vcf_to_msa_args="${SGE_ARGS} -hold_jid $jobid_calling"
	run_vcf_to_tsv=$( qsub $vcf_to_tsv_args -t 1-$sample_count -N $JOBNAME.VCF_TO_TSV $vcf_to_tsv_cmd)
        jobid_vcf_to_tsv=$( echo $run_vcf_to_tsv | cut -d ' ' -f3 | cut -d '.' -f1 )
	
	tsv_to_msa_args="${SGE_ARGS} -hold_jid $jobid_vcf_to_tsv"
        run_tsv_to_msa=$( qsub $tsv_to_msa_args -t 1-$sample_count -N $JOBNAME.TSV_TO_MSA $tsv_to_msa_cmd)
        jobid_tsv_to_msa=$( echo $run_vcf_to_tsv | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "vcf to msa:$jobid_vcf_to_tsv - $jobid_tsv_to_msa\n" >> $output_dir/logs/jobids.txt

	else

	echo "Running variant calling on sample"
	run_vcf_to_tsv=$($vcf_to_tsv_cmd)
	run_tsv_to_msa=$($tsv_to_msa_cmd)
	fi
fi

