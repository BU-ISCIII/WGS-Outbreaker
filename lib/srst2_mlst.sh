#!/bin/bash
## Author: A.Hernandez
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to run srst2 for mlst\n"
        echo -e "Usage: srst2_mlst.sh input_dir output_dir srst2_mlst_fasta srst2_mlst_profile samples_list FASTQ_R1_list FASTQ_R2_list mlst_list srst2_delim"
        exit
fi

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
srst2_delim=$9

if [ "$use_sge" = "1" ]; then
        sample_count=$SGE_TASK_ID
else
        sample_count=${10}
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R1=$( echo $fastq_files_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R2=$( echo $fastq_files_R2_list | tr ":" "\n" | head -$sample_count | tail -1)
mlst_name=$( echo $mlst_list | tr ":" "\n" | head -$sample_count | tail -1)

echo -e "Running mlst.sh for $sample"

srst2 --input_pe $dir/$sample/$fastq_files_R1 $dir/$sample/$fastq_files_R2 --forward trimmed_R1 --reverse trimmed_R2 --log --output $output_dir/$sample/$mlst_name --mlst_db $srst2_db_path_mlst_db --mlst_definitions $srst2_db_path_mlst_definitions --mlst_delimiter "$srst2_delim"

echo -e "mlst.sh for $sample finished"
