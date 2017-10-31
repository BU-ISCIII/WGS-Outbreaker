#!/bin/bash
## Author A.Hernandez
## Usage: run_identification_ST.sh ....

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

CONFIG_FILE=$1


#Execure processing_config.sh
source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

##create directories
mkdir -p $output_dir/kmerfinder

jobid_compress=$( cat $output_dir/logs/jobids.txt | grep -w "COMPRESS_FILE" | cut -d ':' -f2 )

if [ -d $output_dir/QC/trimmomatic ]; then
	dir=$output_dir/QC/trimmomatic
	if [ "$use_sge" = "1" ]; then
		concat_arg="${SGE_ARGS} -pe openmp $threads -hold_jid ${jobid_compress}"
	fi
else
	dir=$output_dir/raw
fi


if [ $trimming == "YES" ]; then

	concat_cmd="$SCRIPTS_DIR/concat.sh \
	$dir \
	$output_dir/kmerfinder \
	$samples \
	$compress_paired_R1_list \
        $compress_paired_R2_list \
	$concatFastq_list"

else
	concat_cmd="$SCRIPTS_DIR/concat.sh \
        $dir \
        $output_dir/kmerfinder \
        $samples \
	$fastq_R1_list \
	$fastq_R2_list \
        $concatFastq_list"
fi

kmerfinder_cmd="$SCRIPTS_DIR/kmerfinder.sh \
	$output_dir/kmerfinder \
	$output_dir/kmerfinder \
	$samples \
	$concatFastq_list \
        $kmerfinderST_list \
	$kmerfinder_path \
	$bact_db_path"

if [ $kmerfinder == "YES" ]; then
	if [ "$use_sge" = "1" ]; then
		concat_files=$( qsub $concat_arg -t 1-$sample_count -N $JOBNAME.CONCATFILES $concat_cmd)
		jobid_concat=$( echo $concat_files | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "CONCAT_FILES:$jobid_concat\n" >> $output_dir/logs/jobids.txt 
			
		kmerfinder_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_concat"
		kmerfinder_qsub=$( qsub $kmerfinder_arg -t 1-$sample_count -N $JOBNAME.KMERFINDER $kmerfinder_cmd)
		jobid_kmerfinder=$( echo $kmerfinder_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "KMERFINDER:$jobid_kmerfinder\n" >> $output_dir/logs/jobids.txt
	else
		for count in  `seq 1 $sample_count`; do
			echo "Running concat files on sample $count"
			concat_files=$($concat_cmd $count)
			echo "Running kmerfinder on sample $count"
			kmerfinder_run=$($kmerfinder_cmd $count)
		done
	fi	
fi



