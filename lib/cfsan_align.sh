#!/bin/bash
##Author:A.Hernandez
## Usage : alignmen_cfsan.sh ...


# Test whether the script is being executed with sge or not.
if [ -z $SGE_TASK_ID ]; then
        use_sge=0
else
        use_sge=1
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
#set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

echo `date`

#VARIABLES

dir=$1
output_dir=$2
samples=$3
fastq_files_R1_list=$4
fastq_files_R2_list=$5
cfsan_ref_path=$6

if [ "$use_sge" = "1" ]; then
        sample_count=$SGE_TASK_ID 
else
        sample_count=$7
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R1=$( echo $fastq_files_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R2=$( echo $fastq_files_R2_list | tr ":" "\n" | head -$sample_count | tail -1)


echo -e "Running aligment cfsan for $sample \n"

cfsan_snp_pipeline map_reads -f $cfsan_ref_path $output_dir/samples/$sample/$fastq_files_R1 $output_dir/samples/$sample/$fastq_files_R2
