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

CONFIG_FILE=$1

#Execure processing_config.sh
source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

#Folder creation

mkdir -p $output_dir/CFSAN/samples
echo "Directory for CFSAN created"

jobid_compress=$(cat $output_dir/logs/jobids.txt | grep -w "COMPRESS_FILE" | cut -d ':' -f2 )

if [ -d $output_dir/QC/trimmomatic ]; then
	dir=$output_dir/QC/trimmomatic
	if [ "$use_sge" = "1" ]; then
		cfsan_arg="${SGE_ARGS} -pe openmp $threads -hold_jid ${jobid_compress}"
	fi
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


# Creation sampleDirectories.txt file

sampledirectories_cmd="$SCRIPTS_DIR/cfsan_create_sample_directories.sh \
	$samples \
	$output_dir/CFSAN"

for count in `seq 1 $sample_count`; do
	echo "creation sampleDirectories.txt file"
        run_sample_directories=$($sampledirectories_cmd $count)

done


create_ln_cmd="$SCRIPTS_DIR/cfsan_create_ln.sh \
	$dir \
	$output_dir \
	$samples \
	$trimming \
	$fastq_files_R1_list \
	$fastq_files_R2_list"

aling_sample_to_reference_cmd="$SCRIPTS_DIR/cfsan_align.sh \
	$dir \
	$output_dir/CFSAN \
	$samples \
	$fastq_files_R1_list \
        $fastq_files_R2_list \
	$ref_path"

picard_sort_sam_cmd="$SCRIPTS_DIR/cfsan_sort_sam.sh \
	$threads \
	$output_dir/CFSAN/samples \
	$samples \
	$sort_sam_list \
	$picard_path"

picard_mark_duplicate_cmd="$SCRIPTS_DIR/picard_duplicates.sh \
	$output_dir/CFSAN/samples \
	$samples
	$sort_sam_list \
	$dedup_sam_list \
	$picard_path"

gatk_add_or_replace_group_cmd="$SCRIPTS_DIR/cfsan_readgroup.sh \
	$threads
	$output_dir/CFSAN/samples \
	$samples
	$dedup_sam_list \
	$dedup_bam_list
	$platform
	$picard_path"


samtool_index_cmd="$SCRIPTS_DIR/cfsan_samtool_index.sh \
	$threads \
	$output_dir/CFSAN/samples \
	$samples \
	$dedup_bam_list"

gatk_realinger_target_creator_cmd="$SCRIPTS_DIR/cfsan_target_creator.sh \
	$threads
	$output_dir/CFSAN/samples \
	$samples
	$dedup_bam_list
	$intervals_list
	$gatk_path
	$ref_path"

gatk_indel_realigner_cmd="$SCRIPTS_DIR/cfsan_indel_realigner.sh \
	$threads \
	$output_dir/CFSAN/samples \
	$samples \
	$intervals_list \
	$dedup_bam_list \
	$unsorted_bam_list \
	$gatk_path \
	$ref_path"

cfsan_call_sites_cmd="$SCRIPTS_DIR/cfsan_call_sites.sh \
	$output_dir/CFSAN/samples \
	$samples
	$unsorted_bam_list
	$cfsan_ref_path"

cfsan_snp_filter_cmd="$SCRIPTS_DIR/cfsan_filter_regions.sh \
	$output_dir/CFSAN \
	$cfsan_ref_path"


cfsan_snplist_cmd="$SCRIPTS_DIR/cfsan_snp_list.sh \
	$output_dir/CFSAN"

cfsan_call_consensus_cmd="$SCRIPTS_DIR/cfsan_call_consensus.sh \
	$output_dir/CFSAN \
	$samples"

cfsan_create_snpmatrix_cmd="$SCRIPTS_DIR/cfsan_snp_matrix.sh \
	$output_dir/CFSAN"

cfsan_snp_reference_cmd="$SCRIPTS_DIR/cfsan_snp_reference.sh \
	$output_dir/CFSAN \
	$cfsan_ref_path"

cfsan_collect_metrics_cmd="$SCRIPTS_DIR/cfsan_metrics.sh \
	$output_dir/CFSAN \
	$samples \
	$cfsan_ref_path"

cfsan_combine_metrics_cmd="$SCRIPTS_DIR/cfsan_combine_metrics.sh \
	$output_dir/CFSAN" 
	
cfsan_merge_vcf_cmd="$SCRIPTS_DIR/cfsan_merge_vcf.sh \
        $output_dir/CFSAN"


cfsan_matrix_distance_cmd="$SCRIPTS_DIR/cfsan_snp_distance.sh \
	$output_dir/CFSAN"

cfsan_clean_cmd=


if [ $cfsan == "YES" ]; then
	if [ "$use_sge" = "1" ]; then
	
	createln_cfsan_qsub=$( qsub $cfsan_arg -t 1-$sample_count -N $JOBNAME.CREATELN $create_ln_cmd)
        jobid_ln_cfsan=$( echo $createln_cfsan_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SYM_LN_FILES:$jobid_ln_cfsan\n" >> $output_dir/logs/jobids.txt
	
	cfsan_align_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_ln_cfsan"
	align_cfsan_qsub=$( qsub $cfsan_align_arg -t 1-$sample_count -N $JOBNAME.ALIGN $aling_sample_to_reference_cmd)
	jobid_align_cfsan=$( echo $align_cfsan_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "ALIGN_FILES:$jobid_align_cfsan\n" >> $output_dir/logs/jobids.txt


	cfsan_callsites_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_align_cfsan"
	cfsan_call_sites_qsub=$( qsub $cfsan_callsites_arg -t 1-$sample_count -N $JOBNAME.CALL.SITES $cfsan_call_sites_cmd)
	jobid_cfsan_call_sites=$( echo $cfsan_call_sites_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "CALL_SITES_FILES:$jobid_cfsan_call_sites\n" >> $output_dir/logs/jobids.txt

	cfsan_snpfil_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_call_sites"
	cfsan_snp_filter_qsub=$( qsub $cfsan_snpfil_arg -N $JOBNAME.SNP.FILTER $cfsan_snp_filter_cmd)
        jobid_cfsan_snpfil=$( echo $cfsan_snp_filter_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SNP_FILTER_FILES:$jobid_cfsan_snpfil\n" >> $output_dir/logs/jobids.txt

	cfsan_snplist_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_snpfil"
	cfsan_snplist_qsub=$( qsub $cfsan_snplist_arg -N $JOBNAME.SNP.LIST $cfsan_snplist_cmd)
        jobid_cfsan_snplist=$( echo $cfsan_snplist_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SNP_LIST_FILES:$jobid_cfsan_snplist\n" >> $output_dir/logs/jobids.txt

	cfsan_callconsen_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_snplist"
	cfsan_call_consensus_qsub=$( qsub $cfsan_callconsen_arg -t 1-$sample_count -N $JOBNAME.CALL.CONSENSUS $cfsan_call_consensus_cmd)
	jobid_cfsan_call_consensus=$( echo $cfsan_call_consensus_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "CALL_CONSENSUS_FILES:$jobid_cfsan_call_consensus\n" >> $output_dir/logs/jobids.txt

	cfsan_snpma_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_call_consensus"
	cfsan_create_snpmatrix_qsub=$( qsub $cfsan_snpma_arg -N $JOBNAME.SNP.MATRIX $cfsan_create_snpmatrix_cmd)
        jobid_cfsan_create_snpmatrix=$( echo $cfsan_create_snpmatrix_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SNP_MATRIX_FILES:$jobid_cfsan_create_snpmatrix\n" >> $output_dir/logs/jobids.txt

	cfsan_snpref_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_create_snpmatrix"
	cfsan_snp_reference_qsub=$( qsub $cfsan_snpref_arg -N $JOBNAME.SNP.REFERENCE $cfsan_snp_reference_cmd)
        jobid_cfsan_snp_reference=$( echo $cfsan_snp_reference_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SNP_REFERENCE_FILES:$jobid_cfsan_snp_reference\n" >> $output_dir/logs/jobids.txt

	cfsan_collecmet_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_snp_reference"
	cfsan_collect_metrics_qsub=$( qsub $cfsan_collecmet_arg -t 1-$sample_count -N $JOBNAME.COLLECT.METRICS $cfsan_collect_metrics_cmd)
        jobid_cfsan_collect_metrics=$( echo $cfsan_collect_metrics_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "COLLECT_METRICS_FILES:$jobid_cfsan_collect_metrics\n" >> $output_dir/logs/jobids.txt

	cfsan_combinemet_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_collect_metrics"
	cfsan_combine_metrics_qsub=$( qsub $cfsan_combinemet_arg  -N $JOBNAME.COMBINE.METRICS $cfsan_combine_metrics_cmd)
        jobid_cfsan_combine_metrics=$( echo $cfsan_combine_metrics_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "COMBINE_METRICS_FILES:$jobid_cfsan_combine_metrics\n" >> $output_dir/logs/jobids.txt

	cfsan_merge_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_combine_metrics"
	cfsan_merge_vcf_qsub=$( qsub $cfsan_merge_arg -N $JOBNAME.MERGE.VCF $cfsan_merge_vcf_cmd)
        jobid_cfsan_merge_vcf=$( echo $cfsan_merge_vcf_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "MERGE_VCF_FILES:$jobid_cfsan_merge_vcf\n" >> $output_dir/logs/jobids.txt

	cfsan_distancema_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_cfsan_merge_vcf"
	cfsan_matrix_distance_qsub=$( qsub $cfsan_distancema_arg -N $JOBNAME.DISTANCEMA.VCF $cfsan_matrix_distance_cmd)
        jobid_cfsan_matrix_distance=$( echo $cfsan_matrix_distance_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "MATRIX_DISTANCE_FILES:$jobid_cfsan_matrix_distance\n" >> $output_dir/logs/jobids.txt


	else
		for count in `seq 1 $sample_count`; do
		echo "Running cfsan for $sample"
		make_ln=$($create_ln_cmd $count)
		
		echo "CFSAN align for $sample"
		run_cfsan_aling=$($aling_sample_to_reference_cmd $count)

		echo "CFSAN call sites for $sample"
		run_cfsan_call_sites=$($cfsan_call_sites_cmd $count)
		
		done

		echo "CFSAN snp filter"
		run_cfsan_snp_filter=$($cfsan_snp_filter_cmd)

		echo "CFSAN snp list"
		run_cfsan_snplist=$($cfsan_snplist_cmd)

		for count in `seq 1 $sample_count`; do

        		echo "CFSAN call consensus for $sample"
        		run_cfsan_call_consensus=$($cfsan_call_consensus_cmd $count)
		done

		echo "CFSAN snp matrix"
		run_cfsan_snpmatrix=$($cfsan_create_snpmatrix_cmd)

		echo "CFSAN snp reference"
		run_cfsan_snp_reference=$($cfsan_snp_reference_cmd)

		for count in `seq 1 $sample_count`; do
        		echo "CFSAN collect metrics"
        		run_cfsan_collect_metrics=$($cfsan_collect_metrics_cmd $count)
		done
		
		echo "CFSAN combine metrics"		
		run_cfsan_combine_metrics=$($cfsan_combine_metrics_cmd)
		
		echo "CFSAN merge vcf"
		run_cfsan_merge=$($cfsan_merge_vcf_cmd)

		echo "CFSAN matrix distance"
		run_matrix_distance=$($cfsan_matrix_distance_cmd)
	
	fi
fi
