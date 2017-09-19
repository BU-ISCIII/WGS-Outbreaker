#!/bin/bash
##Author:A.Hernandez
##Usage: run_srst2.sh ...

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


#Execure processing_config.sh
source $SCRIPTS_DIR/processing_config.sh

#Create directories

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)

for count in $sample; do
	mkdir -p $output_dir/srst2/$sample
	echo "Directory for $sample created"
	done


jobid_trimmomatic=$( cat $output_dir/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )

if [ -d $output_dir/QC/trimmomatic ]; then
        dir=$output_dir/QC/trimmomatic
        srst2_arg="${SGE_ARGS} -hold_jid ${jobid_trimmomatic}"
else
        dir=$output_dir/raw
fi

if [ $trimming == "YES" ];then
	fastq_files_R1_list=$compress_paired_R1_list
	fastq_files_R2_list=$compress_paired_R2_list
else
	fastq_files_R1_list=$fastq_R1_list
        fastq_files_R2_list=$fastq_R2_list
fi

resistance_cmd="$SCRIPTS_DIR/resistance.sh \
	$dir \
	$output_dir/srst2 \
	$samples \
	$fastq_files_R1_list \
	$fastq_files_R2_list \
	$resistance_list \
	$srst2_db_path_argannot"


plasmid_cmd="$SCRIPTS_DIR/plasmid.sh \
	$dir \
	$output_dir/srst2 \
	$samples \
	$fastq_files_R1_list \
	$fastq_files_R2_list \
	$plasmid_list \
	$srst2_db_path_plasmidfinder"


mlst_cmd="$SCRIPTS_DIR/mlst.sh \
	$dir \
	$output_dir/srst2 \
	$srst2_db_path_mlst_db \
	$srst2_db_path_mlst_definitions \
	$samples \
	$fastq_files_R1_list \
        $fastq_files_R2_list \
	$mlst_list"


summary_cmd="$SCRIPTS_DIR/summary_srst2.sh \
$output_dir/srst2 
$output_dir/srst2" 

if [ $srst2 == "YES" ];then
	if [ "$use_sge" = "1" ]; then
		resistance_qsub=$( qsub $srst2_arg -N $JOBNAME.RESISTANCE $resistance_cmd )
		jobid_resistance=$( echo $resistance_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "RESISTANCE FILES:$jobid_resistance\n" >> $output_dir/logs/jobids.txt

		plasmids_qsub=$( qsub $srst2_arg -N $JOBNAME.PLASMIDS $plasmid_cmd )
                jobid_plasmid=$( echo $plasmids_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "PLASMID FILES:$jobid_plasmid\n" >> $output_dir/logs/jobids.txt
	
		mlst_qsub=$( qsub $srst2_arg -N $JOBNAME.MLST $mlst_cmd )
                jobid_mlst=$( echo $mlst_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "MLST FILES:$jobid_mlst\n" >> $output_dir/logs/jobids.txt
		
		summary_qsub=$( qsub $srst2_arg -N $JOBNAME.SUMMARY $summary_cmd )
                jobid_mlst=$( echo $summary | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "SUMAMRY FILES:$jobid_summary\n" >> $output_dir/logs/jobids.txt

		
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

