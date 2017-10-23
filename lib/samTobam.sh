#!/bin/bash
## Author S. Monzon
## version v2.0 

# Test whether the script is being executed with sge or not.                                                                                                                                                                        
if [ -z $SGE_TASK_ID ]; then                                                                                                                                                                                                        
  	use_sge=0                                                                                                                                                                                                                        
else                                                                                                                                                                                                                                
  	use_sge=1                                                                                                                                                                                                                        
fi                                                                                                                                                                                                                                                                                                                                                        
## Usage                                                                                                                                                                                                                  
if [ $# != 14 -a "$use_sge" == "1" ]; then                                                                                                                                                                                           
  	echo "usage: ............"                                                                                                                                                                                                       
  	exit                                                                                                                                                                                                                             
elif [ $# != 15 -a "$use_sge" == "0" ]; then                                                                                                                                                                                         
  	echo "usage: ............"                                                                                                                                                                                                       
   	exit                                                                                                                                                                                                                             
fi                                                                                                                                                                                                                                  
# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

echo `date`                                                                                                                                                                                                                              
# VARIABLES
             
input_dir=$1
output_dir=$2
samples=$3
mappingArray_sam_list=$4
mappingArray_bam_list=$5
mappingArray_sorted_list=$6
mappingArray_rg_list=$7
picard_path=$8
platform=$9
model=${10}
date_run=${11}
library=${12}
sequencing_center=${13}
run_platform=${14}



if [ "$use_sge" = "1" ]; then                                                                                                                                                                                                        
  	sample_count=$SGE_TASK_ID                                                                                                                                                                                                      
else                                                                                                                                                                                                                                 
  	sample_count=${15}                                                                                                                                                                                                                 
fi                                                                                                                                                                                                                                                                                                                                                  
sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
mappingArray_sam=$( echo $mappingArray_sam_list | tr ":" "\n" | head -$sample_count | tail -1)
mappingArray_bam=$( echo $mappingArray_bam_list | tr ":" "\n" | head -$sample_count | tail -1)
mappingArray_sorted=$( echo $mappingArray_sorted_list | tr ":" "\n" | head -$sample_count | tail -1)
mappingArray_rg=$( echo $mappingArray_rg_list | tr ":" "\n" | head -$sample_count | tail -1)  

mkdir -p $output_dir/$sample

samtools view -bS $input_dir/$sample/$mappingArray_sam -o $output_dir/$sample/$mappingArray_bam
samtools sort -o $output_dir/$sample/$mappingArray_sorted -T $output_dir/$sample/$mappingArray_sorted $output_dir/$sample/$mappingArray_bam
          
java $JAVA_RAM -jar $picard_path/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT INPUT=$output_dir/$sample/$mappingArray_sorted OUTPUT=$output_dir/$sample/$mappingArray_rg RGID=$date_run-$library-$model-$platform-$sequencing_center RGLB=$library RGPL=$platform RGSM=$sample RGPU=$run_platform RGDT=$date_run RGCN=$sequencing_center

rm -r $input_dir/$sample                                                                                                                             
# Se indexa el fichero BAM con samtools.                                                                                                           
samtools index $output_dir/$sample/$mappingArray_rg
