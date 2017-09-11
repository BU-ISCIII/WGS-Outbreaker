#!/bin/bash
## Author: A.Hernandez
## Usage: resistance.sh


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


# VARIABLES

dir=$1
output_dir=$2
samples=$3
fastq_files_R1_list=$4
fastq_files_R2_list=$5
resistance_list=$6
srst2_db_path_argannot=$7

if [ "$use_sge" = "1" ]; then
        sample_count=$sge_task_id
else
        sample_count=$8
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R1=$( echo $fastq_files_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R2=$( echo $fastq_files_R2_list | tr ":" "\n" | head -$sample_count | tail -1)
resistance_name=$( echo $resistance_list | tr ":" "\n" | head -$sample_count | tail -1)

echo -e "Running resistance.sh for $sample \n\n"

srst2 --input_pe $dir/$sample/$fastq_files_R1 $dir/$sample/$fastq_files_R2 --forward trimmed_R1 --reverse trimmed_R2 --log --output $output_dir/$sample/$resistance_name --gene_db $srst2_db_path_argannot

echo -e "Resistance.sh for $sample finished \n\n"
