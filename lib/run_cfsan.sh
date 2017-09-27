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


# Creation sampleDirectories.txt file

sampledirectories_cmd="$SCRIPTS_DIR/create_sample_directories.sh \
	$samples \
	$output_dir/CFSAN"

for count in `seq 1 $sample_count`; do
	echo "creation sampleDirectories.txt file"
        run_sample_directories=$($sampledirectories_cmd $count)

done


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
	$threads
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
	align_cfsan_qsub=$( qsub $cfsan_arg -N $JOBNAME.ALIGN $aling_sample_to_reference_cmd)
	jobid_align_cfsan=$( echo $align_cfsan_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "ALIGN FILES:$jobid_alig_cfsan\n" >> output_dir/logs/jobids.txt
	
	cfsan_sort_sam_qsub=$( qsub $cfsan_arg -N $JOBNAME.SORT_SAM $picard_sort_sam_cmd)
        jobid_cfsan_sort_sam=$( echo $cfsan_sort_sam_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SORT SAM FILES:$jobid_cfsan_sort_sam\n" >> output_dir/logs/jobids.txt
	
	cfsan_duplicates_qsub=$( qsub $cfsan_arg -N $JOBNAME.SORT_SAM $picard_mark_duplicate_cmd)
        jobid_cfsan_duplicates=$( echo $cfsan_duplicates_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "MARK DUPLICATES FILES:$jobid_cfsan_duplicates\n" >> output_dir/logs/jobids.txt

	cfsan_read_group_qsub=$( qsub $cfsan_arg -N $JOBNAME.READGROUP $gatk_add_or_replace_group_cmd)
	jobid_cfsan_read_group=$( echo $cfsan_read_group_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "READ GROUP FILES:$jobid_cfsan_read_group\n" >> output_dir/logs/jobids.txt
	
	cfsan_samtool_index_qsub=$( qsub $cfsan_arg -N $JOBNAME.CFSAN.INDEX $samtool_index_cmd)
	jobid_cfsan_samtool_index=$( echo $cfsan_samtool_index_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SAMTOOL INDEX FILES:$jobid_cfsan_samtool_index\n" >> output_dir/logs/jobids.txt

	cfsan_target_creator_qsub=$( qsub $cfsan_arg -N $JOBNAME.CFSAN.TARGET.CREATOR $gatk_realinger_target_creator_cmd)
	jobid_cfsan_target_creator=$( echo $cfsan_target_creator_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "TARGET CREATOR FILES:$jobid_cfsan_target_creator\n" >> output_dir/logs/jobids.txt

	cfsan_indel_realigner_qsub=$( qsub $cfsan_arg -N $JOBNAME.INDEL.REALIGNER  $gatk_indel_realigner_cmd)
	jobid_cfsan_indel_realigner=$( echo $cfsan_indel_realigner_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "INDEL REALIGNER FILES:$jobid_cfsan_indel_realigner\n" >> output_dir/logs/jobids.txt

	cfsan_call_sites_qsub=$( qsub $cfsan_arg -N $JOBNAME.CALL.SITES $cfsan_call_sites_cmd)
	jobid_cfsan_call_sites=$( echo $cfsan_call_sites_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "CALL SITES FILES:$jobid_cfsan_call_sites\n" >> output_dir/logs/jobids.txt

	cfsan_snp_filter_qsub=$( qsub $cfsan_arg -N $JOBNAME.SNP.FILTER $cfsan_snp_filter_cmd)
        jobid_cfsan_snp_filter=$( echo $cfsan_snp_filter_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SNP FILTER FILES:$jobid_cfsan_snp_filter\n" >> output_dir/logs/jobids.txt

	cfsan_snplist_qsub=$( qsub $cfsan_arg -N $JOBNAME.SNP.LIST $cfsan_snplist_cmd)
        jobid_cfsan_snplist=$( echo $cfsan_snplist_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SNP LIST FILES:$jobid_cfsan_snplist\n" >> output_dir/logs/jobids.txt

	cfsan_call_consensus_qsub=$( qsub $cfsan_arg -N $JOBNAME.CALL.CONSENSUS $cfsan_call_consensus_cmd)
	jobid_cfsan_call_consensus=$( echo $cfsan_call_consensus_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "CALL CONSENSUS FILES:$jobid_cfsan_call_consensus\n" >> output_dir/logs/jobids.txt

	cfsan_create_snpmatrix_qsub=$( qsub $cfsan_arg -N $JOBNAME.SNP.MATRIX $cfsan_create_snpmatrix_cmd)
        jobid_cfsan_create_snpmatrix=$( echo $cfsan_create_snpmatrix_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SNP MATRIX FILES:$jobid_cfsan_create_snpmatrix\n" >> output_dir/logs/jobids.txt

	cfsan_snp_reference_qsub=$( qsub $cfsan_arg -N $JOBNAME.SNP.REFERENCE $cfsan_snp_reference_cmd)
        jobid_cfsan_snp_reference=$( echo $cfsan_snp_reference_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "SNP REFERENCE FILES:$jobid_cfsan_snp_reference\n" >> output_dir/logs/jobids.txt

	cfsan_collect_metrics_qsub=$( qsub $cfsan_arg -N $JOBNAME.COLLECT.METRICS $cfsan_collect_metrics_cmd)
        jobid_cfsan_collect_metrics=$( echo $cfsan_collect_metrics_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "COLLECT METRICS FILES:$jobid_cfsan_collect_metrics\n" >> output_dir/logs/jobids.txt

	cfsan_combine_metrics_qsub=$( qsub $cfsan_arg -N $JOBNAME.COMBINE.METRICS $cfsan_combine_metrics_cmd)
        jobid_cfsan_combine_metrics=$( echo $cfsan_combine_metrics_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "COMBINE METRICS FILES:$jobid_cfsan_combine_metrics\n" >> output_dir/logs/jobids.txt

	cfsan_merge_vcf_qsub=$( qsub $cfsan_arg -N $JOBNAME.MERGE.VCF $cfsan_merge_vcf_cmd)
        jobid_cfsan_merge_vcf=$( echo $cfsan_merge_vcf_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "MERGE VCF FILES:$jobid_cfsan_merge_vcf\n" >> output_dir/logs/jobids.txt

	cfsan_matrix_distance_qsub=$( qsub $cfsan_arg -N $JOBNAME.MERGE.VCF $cfsan_matrix_distance_cmd)
        jobid_cfsan_matrix_distance=$( echo $cfsan_matrix_distance_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "MATRIX DISTANCE FILES:$jobid_cfsan_matrix_dsitance\n" >> output_dir/logs/jobids.txt



	else
		for count in `seq 1 $sample_count`; do
		echo "Running cfsan for $sample"
		make_ln=$($create_ln_cmd $count)
		
		echo "CFSAN align for $sample"
		run_cfsan_aling=$($aling_sample_to_reference_cmd $count)
		
		echo "CFSAN sort sam for $sample"
		run_cfsan_sort_sam=$($picard_sort_sam_cmd $count)
		
		echo "CFSAN mark duplicates for $sample"
		run_cfsan_duplicates=$($picard_mark_duplicate_cmd $count)

		echo "CFSAN read group for $sample"
                run_cfsan_read_group=$($gatk_add_or_replace_group_cmd $count)
                
		echo "CFSAN samtool index for $sample"
                run_cfsan_samtool_index=$($samtool_index_cmd $count)

		echo "CFSAN target creator for $sample"
		run_cfsan_target_creator=$($gatk_realinger_target_creator_cmd $count)
		
		echo "CFSAN indel realigner for $sample"
		run_cfsan_indel_realigner=$($gatk_indel_realigner_cmd $count)

		echo "CFSAN call sites for $sample"
		run_cfsan_call_sites=$($cfsan_call_sites_cmd $count)
		
		done
	
	fi
fi

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

echo "CFSAN merge VCF"
run_cfsan_merge_vcf=$($cfsan_merge_vcf_cmd)

echo "CFSAN matrix distance"
run_cfsan_matrix_distance=$($cfsan_matrix_distance_cmd)

