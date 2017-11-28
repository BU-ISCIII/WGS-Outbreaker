#!/bin/bash
## Author:A.Hernandez
## version v2.0

if [ $# -eq 0 ];then
	echo -e "\nScript to create symbolic links to Fastq files for cfsan pipeline\n"
	echo -e "Usage: cfsan_create_ln.sh input_dir output_dir samples_list trimming_step FastqR1_list FastqR2_list"
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

#VARIABLES

dir=$1
output_dir=$2
samples=$3
trimming=$4
fastq_files_R1_list=$5
fastq_files_R2_list=$6

if [ "$use_sge" = "1" ]; then
        sample_count=$SGE_TASK_ID
else
        sample_count=$7
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R1=$( echo $fastq_files_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
fastq_files_R2=$( echo $fastq_files_R2_list | tr ":" "\n" | head -$sample_count | tail -1)

# Folder creation

mkdir -p $output_dir/CFSAN/samples/$sample
echo "creation  CFSAN folder for $sample"

if [ $trimming == "YES" ]; then
	for count in $samples; do
	ln -fs $dir/$sample/$fastq_files_R1 $output_dir/CFSAN/samples/$sample/$fastq_files_R1
	ln -fs $dir/$sample/$fastq_files_R2 $output_dir/CFSAN/samples/$sample/$fastq_files_R2
	done
else
	for count in $samples; do
	ln -fs $dir/$sample/$fastq_files_R1 $output_dir/CFSAN/samples/$sample/$fastq_files_R1
	ln -fs $dir/$sample/$fastq_files_R2 $output_dir/CFSAN/samples/$sample/$fastq_files_R2
	done
fi
