#!/bin/bash
## Author S. Monzon                                                                                                  ## version v2.0                                                                                                                                                                                                                          
                                                                                                                     # Test whether the script is being executed with sge or not.                                                                                                                                                                       
if [ -z $SGE_TASK_ID ]; then                                                                                         	use_sge=0                                                                                                                                                                                                            
else                                                                                                                    use_sge=1                                                                                                                                                                                                                     
fi                                                                                                                                                                                                                                                                                                                                           
## Usage                                                                                                                                                                                                                          
if [ $# != 5 -a "$use_sge" == "1" ]; then                                                                                                                                                                                        
   	echo "usage: ............"                                                                                                                                                                 
   	exit                                                                                                                                                                                                                     
elif [ $# != 6 -a "$use_sge" == "0" ]; then                                                                                                                                                                                   
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
# Variables
input_dir=$1
samples=$2
mappingArray_rg_list=$3
duplicateBamArray_list=$4
picard_path=$5
                                                                                                                                                                                                                              
if [ "$use_sge" = "1" ]; then                                                                                                                                                                                               
   	sample_count=$sge_task_id                                                                      
else                                                                                                        
   	sample_count=$6                                                                                     
fi
                                                                                                                                                                                              
sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
mappingArray_rg=$( echo $mappingArray_rg_list | tr ":" "\n" | head -$sample_count | tail -1) 
duplicateBamArray=$( echo $duplicateBamArray_list | tr ":" "\n" | head -$sample_count | tail -1) 

echo -e "- Duplicate Filter: $sample"                                                                                                                                                                                                                                                                                       
java $JAVA_RAM -jar $picard_path/picard.jar MarkDuplicates ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT="$input_dir/$sample/$mappingArray_rg=" OUTPUT="$input_dir/$sample/$duplicateBamArray" METRICS_FILE="$input_dir/$sample/$sample.duplicates.stats" TMP_DIR="$TEMP"

samtools index $input_dir/$sample/$duplicateBamArray                                                      
