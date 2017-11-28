#!/bin/bash
## Author S. Monzon
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to run FastQC\n"
        echo -e "Usage: fastqc.sh theads input_dir output_dir samples_list FastqR1_list FastqR2_list"
        exit
fi


#  Test whether the script is being executed with sge or not.
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

## VARIABLES

threads=$1
input_dir=$2
output_dir=$3
samples=$4
fastq_R1_list=$5
fastq_R2_list=$6

if [ "$use_sge" = "1" ]; then
	sample_count=$SGE_TASK_ID
else
	sample_count=$7
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
fastq_R1=$( echo $fastq_R1_list | tr ":" "\n" | head -$sample_count | tail -1) 
fastq_R2=$( echo $fastq_R2_list | tr ":" "\n" | head -$sample_count | tail -1)  

echo -e "Running FastQC for $sample....\n"

# Results folder per sample creation
mkdir -p $output_dir/QC/fastqc/$sample

fastqc --noextract -o $output_dir/QC/fastqc/$sample -t $threads $input_dir/$sample/$fastq_R1 $input_dir/$sample/$fastq_R2	 

echo -e "Preprocessing for $sample finished \n\n"
