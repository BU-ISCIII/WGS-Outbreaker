#!/bin/bash
## Author:A.Hernandez
## version v2.0

if [ $# -eq 0  ]; then
	echo -e "\nExecute SRST2 for resisitance, plasmids and mlst\n"
        echo "Usage: run_srst2.sh <config.file>"
        exit
fi

#Execute processing_config.sh

CONFIG_FILE=$1

# Check if run_outbreak_wgs.sh was execute
if [ -z $SCRIPTS_DIR ]; then
        SCRIPTS_DIR=$( cat $CONFIG_FILE | grep -w 'SCRIPTS_DIR' | cut -d '=' -f2 )
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

# Or other runner was execute
else
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"
fi


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


#Folder creation

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)

for count in $sample; do
	mkdir -p $output_dir/srst2/$sample
	echo "Directory for $sample created"
	done

# Get jobid for compress files step
jobid_compress=$( cat $output_dir/logs/jobids.txt | grep -w "COMPRESS_FILE" | cut -d ':' -f2 )

# Check if trimmomatic was executed for input files
if [ -d $output_dir/QC/trimmomatic ]; then
        dir=$output_dir/QC/trimmomatic
	if [ "$use_sge" = "1" ]; then
        	srst2_arg="${SGE_ARGS} -pe openmp $threads -hold_jid ${jobid_compress}"
	fi
else
        dir=$output_dir/raw
fi

# Set fastq files names
if [ $trimming == "YES" ];then
	fastq_files_R1_list=$compress_paired_R1_list
	fastq_files_R2_list=$compress_paired_R2_list
else
	fastq_files_R1_list=$fastq_R1_list
        fastq_files_R2_list=$fastq_R2_list
fi

resistance_cmd="$SCRIPTS_DIR/srst2_resistance.sh \
	$dir \
	$output_dir/srst2 \
	$samples \
	$fastq_files_R1_list \
	$fastq_files_R2_list \
	$resistance_list \
	$srst2_db_path_argannot"

plasmid_cmd="$SCRIPTS_DIR/srst2_plasmid.sh \
	$dir \
	$output_dir/srst2 \
	$samples \
	$fastq_files_R1_list \
	$fastq_files_R2_list \
	$plasmid_list \
	$srst2_db_path_plasmidfinder"

mlst_cmd="$SCRIPTS_DIR/srst2_mlst.sh \
	$dir \
	$output_dir/srst2 \
	$srst2_db_path_mlst_db \
	$srst2_db_path_mlst_definitions \
	$samples \
	$fastq_files_R1_list \
        $fastq_files_R2_list \
	$mlst_list \
	$srst2_delim"

summary_cmd="$SCRIPTS_DIR/srst2_summary.sh \
	$output_dir/srst2 \
	$output_dir/srst2" 

# Execute SRST2
if [ $srst2 == "YES" ];then
	#In HPC
	if [ "$use_sge" = "1" ]; then
		resistance_qsub=$( qsub $srst2_arg -t 1-$sample_count -N $JOBNAME.RESISTANCE $resistance_cmd )
		jobid_resistance=$( echo $resistance_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "RESISTANCE_FILES:$jobid_resistance\n" >> $output_dir/logs/jobids.txt

		plasmids_qsub=$( qsub $srst2_arg -t 1-$sample_count -N $JOBNAME.PLASMIDS $plasmid_cmd )
                jobid_plasmid=$( echo $plasmids_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "PLASMID_FILES:$jobid_plasmid\n" >> $output_dir/logs/jobids.txt
	
		mlst_qsub=$( qsub $srst2_arg -t 1-$sample_count -N $JOBNAME.MLST $mlst_cmd )
                jobid_mlst=$( echo $mlst_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "MLST_FILES:$jobid_mlst\n" >> $output_dir/logs/jobids.txt
		
		summary_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_mlst"
		summary_qsub=$( qsub $summary_args -N $JOBNAME.SUMMARY $summary_cmd )
                jobid_summary=$( echo $summary_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "SUMAMRY_FILES:$jobid_summary\n" >> $output_dir/logs/jobids.txt

	#Or local	
	else
		for count in `seq 1 $sample_count`; do
			echo "Running resistance on sample $count"
			resistance_run=$($resistance_cmd $count)
			echo "Running plasmid on sample $count"
               		plasmid_run=$($plasmid_cmd $count)
			echo "Running mlst on sample $count"
        		mlst_run=$($mlst_cmd $count)
		done
			
		echo "Running summary"
		summary_run=$($summary_cmd)
		
	fi
fi

