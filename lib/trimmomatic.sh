#!/bin/bash                                                      
## Author: S.Monzon                       
# usage: trimmomatic.sh .....
                                              

####################                                             
###   Commands  ####                                             
####################         


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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
## Usage                                                                                                                                                                                                                                 
echo `date`                                                                                                                                                                                                                        
                                                                                                                                                                                                                                    
## VARIABLES
                                                                                                 
threads=$1
input_dir=$2
output_dir=$3
samples=$4
fastq_R1_list=$5
fastq_R2_list=$6
trimmedFastqArray_paired_R1_list=$7
trimmedFastqArray_paired_R2_list=$8
trimmedFastqArray_unpaired_R1_list=$9
trimmedFastqArray_unpaired_R2_list=${10}
trim_args=${11}
trimmomatic_version=${12}
trimmomatic_path=${13}

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                                      
 	sample_count=$SGE_TASK_ID                                                                                                                                                                                               
else                                                                                                                                                                                                                               
 	sample_count=${14}                                                                                                                                                                                                               
fi                                                                                                                                                                                                                                                                                                                            
sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
fastq_R1=$( echo $fastq_R1_list | tr ":" "\n" | head -$sample_count | tail -1) 
fastq_R2=$( echo $fastq_R2_list | tr ":" "\n" | head -$sample_count | tail -1)
trimmedFastqArray_paired_R1=$( echo $trimmedFastqArray_paired_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
trimmedFastqArray_paired_R2=$( echo $trimmedFastqArray_paired_R2_list | tr ":" "\n" | head -$sample_count | tail -1)
trimmedFastqArray_unpaired_R1=$( echo $trimmedFastqArray_unpaired_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
trimmedFastqArray_unpaired_R2=$( echo $trimmedFastqArray_unpaired_R2_list | tr ":" "\n" | head -$sample_count | tail -1)   
trim_args=$(echo $trim_args | tr "_" " ")

echo -e "Running Trimmomatic for $sample....\n"                                                                                                                                                                                                                            
# Results folder per sample creation                                                                                                                                                                                               
mkdir -p $output_dir/QC/trimmomatic/$sample                                                                                                                                                                                         

java -jar $trimmomatic_path/trimmomatic-$trimmomatic_version.jar PE -phred33 $input_dir/$fastq_R1 $input_dir/$fastq_R2 $output_dir/QC/trimmomatic/$sample/$trimmedFastqArray_paired_R1 $output_dir/QC/trimmomatic/$sample/$trimmedFastqArray_unpaired_R1 $output_dir/QC/trimmomatic/$sample/$trimmedFastqArray_paired_R2 $output_dir/QC/trimmomatic/$sample/$trimmedFastqArray_unpaired_R2 $trim_args
	                                                                                                             echo -e "Trimmomatic for $sample finished \n\n"
