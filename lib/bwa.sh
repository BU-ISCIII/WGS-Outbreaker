#!/bin/bash
## Author A. Hernandez
## version v2.0                                                                                                                                                                                                                                   
# Test whether the script is being executed with sge or not.
if [ -z $sge_task_id ]; then                                                                                                                                                                                                       
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
if [ $# != 8 -a "$use_sge" == "1" ]; then                                                                                                                                                                                          
 	echo "usage: ............"                                                                                      
 	exit                                                                                                                                                                                                                           
elif [ $# != 9 -a "$use_sge" == "0" ]; then                                                                                                                                                                                        
 	echo "usage: ............"                                                                             
  	exit                                                                                                                                                                                                                           
fi                                                                                                                                                                                                                             
echo `date` 

# VARIABLES

threads=$1
input_dir=$2
output_dir=$3
samples=$4
fastq_R1_list=$5
fastq_R2_list=$6
mappingArray_sam_list=$7
ref_path=$8

if [ "$use_sge" = "1" ]; then                                                      
 	sample_count=$sge_task_id                                                  
else                                                                               
 	sample_count=${9}                                                               
fi                                                                                 
                                                                                                                                                                        
sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)       
fastq_R1=$( echo $fastq_R1_list | tr ":" "\n" | head -$sample_count | tail -1)   
fastq_R2=$( echo $fastq_R2_list | tr ":" "\n" | head -$sample_count | tail -1)   
mappingArray_sam=$( echo $mappingArray_sam_list | tr ":" "\n" | head -$sample_count | tail -1)    

mkdir -p $output_dir/$sample

echo -e "2) Alignment using BWA: $sample"                                                                                                                                                                                                                                                                             
echo "Alignment PE"
# Use bwa mem with pair end data align 70bp-1Mbp query sequences

bwa mem -t $threads $ref_path $input_dir/$sample/$fastq_R1 $input_dir/$sample/$fastq_R2 > $output_dir/$sample/$mappingArray_sam     
