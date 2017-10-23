#!/bin/bash
## Author: A.Hernandez
## Usage: mlst.sh ..

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


# VARAIBLES

dir=$1
output_dir=$2
srst2_db_path_mlst_db=$3
srst2_db_path_mlst_definitions=$4
samples=$5
fastq_files_R1_list=$6
fastq_files_R2_list=$7
mlst_list=$8

if [ "$use_sge" = "1" ]; then
        sample_count=$SGE_TASK_ID
else
        sample_count=$9
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R1=$( echo $fastq_files_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R2=$( echo $fastq_files_R2_list | tr ":" "\n" | head -$sample_count | tail -1)
mlst_name=$( echo $mlst_list | tr ":" "\n" | head -$sample_count | tail -1)

echo -e "Running mlst.sh for $sample"

srst2 --input_pe $dir/$sample/$fastq_files_R1 $dir/$sample/$fastq_files_R2 --forward trimmed_R1 --reverse trimmed_R2 --log --output $output_dir/$sample/$mlst_name --mlst_db $srst2_db_path_mlst_db --mlst_definitions $srst2_db_path_mlst_definitions --mlst_delimiter '-'

echo -e "mlst.sh for $sample finished"

