#!/bin/bash
## Author: A. Hernandez
# Help 
# usage: run_cfsan.sh ....
#

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#Execure processing_config.sh
source $SCRIPTS_DIR/processing_config.sh

#Folder creation

mkdir -p $output_dir/CFSAN/samples
echo "Directory for CFSAN created"

jobid_trimmomatic=$(cat $output_dir/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )

if [ -d $output_dir/QC/trimmomatic ]; then
	dir=$output_dir/QC/trimmomatic
	cfsan_args="${SGE_ARGS} -pe orte $threads -hold_jid ${jobid_trimmomatic}"
else
  	dir=$output_dir/raw
fi

if [ $trimming == "YES" ]; then
	fastq_files_R1_list=$compress_paired_R1_list
	fastq_files_R2_list=$compress_paired_R2_list
else
	fastq_files_R1_list=$fastq_R1_list
        fastq_files_R2_list=$fastq_R2_list
fi

create_ln_cmd="$SCRIPTS_DIR/create_ln.sh \
	$dir \
	$output_dir \
	$samples \
	$trimming \
	$fastq_files_R1_list \
	$fastq_files_R2_list"

prepare_reference_cmd=

aling_sample_to_reference_cmd="$SCRIPTS_DIR/alingnmen_cfsan.sh \
	$threads \
	$dir \
	$output_dir/CFSAN \
	$samples \
	$fastq_files_R1_list \
        $fastq_files_R2_list \
	$cfsan_ref_path"

picard_sort_sam_cmd=

picard_mark_duplicate_cmd=

gatk_add_or_replace_group_cmd=

samtool_inex_cmd=

gatk_realinger_target_creator_cmd=

gatk_indel_realigner_cmd=

cfsan_prepsamples_cmd=

cfsan_snp_filer_cmd=

cfsan_snp_filter_fil_cmd=

cfsan_snplist_cmd=

cfsan_snplist_fil_cmd=

cfsan_call_consensus_cmd=

cfsan_call_consensus_fil_cmd=

cfsan_create_snpmatrix_cmd=

sfsan_create_snpmatrix_fil_cmd=

cfsan_snp_reference_cmd=

cfsan_snp_reference_fil_cmd=

cfsan_collect_metrics_cmd=

cfsan_combine_metrics_cmd=

cfsan_matrix_distance_cmd=

cfsan_matrix_distance_fil_cmd=

cfsan_merge_vcf_cmd=

cfsan_clean_cmd=


if [ $cfsan == "YES" ]; then
	if [ "$use_sge" = "1" ]; then
	align_cfsan_qsub=$( qsub $cfsan_arg -N $JOBNAME.ALIGN $aling_sample_to_reference_cmd)
	jobid_align_cfsan=$( echo $align_cfsan_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "ALIGN FILES:$jobid_alig_cfsan\n" >> output_dir/logs/jobids.txt
	
	else
		for count in `seq 1 $sample_count`; do
		echo "Running cfsan for $sample"
		make_ln=$($create_ln_cmd $count)
		
		echo "CFSAN align for $sample"
		run_cfsan_aling=$($aling_sample_to_reference_cmd $count)
		done
	fi
fi
	
