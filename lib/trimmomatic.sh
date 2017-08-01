#!/bin/bash                                                      
                       
# usage: trimmomatic.sh .....
#                                                                

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
                                                                                                                                                                                                                                    
## Usage                                                                                                                                                                                                                           
echo $#

if [ $# != 13 -a "$use_sge" == "1" ]; then                                                                                                                                                                                          
	echo "usage: <input dir> <output dir> <sample names (s1:s2:sn)> <fastq files R1 (f1:f2:fn)> <fastq files R2 (f1:f2:fn)> <trimmed output files paired R1> <trimmed output files paired R2> <trimmed output files unpaired R1> <trimmed output files unpaired R2> <threads> <trim_args> <trimmomatic_version> <trimmomatic path>"                                                                                             
 	exit                                                                                                                                                                                                                           
elif [ $# != 14 -a "$use_sge" == "0" ]; then                                                                                                                                                                                        
 	echo "usage: <input dir> <output dir> <sample names (s1:s2:sn)> <fastq files R1 (f1:f2:fn)> <fastq files R2 (f1:f2:fn)> <trimmed output files paired R1> <trimmed output files paired R2> <trimmed output files unpaired R1> <trimmed output files unpaired R2> <threads> <trim_args> <trimmomatic_version> <trimmomatic path> <sample_number>"                                                                             
  	exit                                                                                                                                                                                                                           
fi                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                    
#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed                     
set -x                                                                                                                                                                                                                             
echo `date`                                                                                                                                                                                                                        
                                                                                                                                                                                                                                    
## Variables                                                                                                                                                                                                                       
INPUT_DIR=$1                                                                                                                                                                                                                       
OUTPUT_DIR=$2                                                                                                                                                                                                                      
THREADS=${10}                                                                                                                                                                                                                         
SAMPLE_NAMES=$3                                                                                                                                                                                                                    
FASTQ_FILES_R1=$4                                                                                                                                                                                                                  
FASTQ_FILES_R2=$5                                                                                                                                                                                                                  
TRIM_ARGS=${11}
TRIM_FILES_PAIRED_R1=$6
TRIM_FILES_PAIRED_R2=$7
TRIM_FILES_UNPAIRED_R1=$8 
TRIM_FILES_UNPAIRED_R2=$9
trimmomatic_version=${12}
TRIMMOMATIC_PATH=${13}

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                                      
 	sample_number=$SGE_TASK_ID                                                                                                                                                                                                     
else                                                                                                                                                                                                                               
 	sample_number=${14}                                                                                                                                                                                                               
fi                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                    
SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                       
FASTQ_R1=$( echo $FASTQ_FILES_R1 | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                   
FASTQ_R2=$( echo $FASTQ_FILES_R2 | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                                                                                                                                                                                                                                                       
trimmedFastqArray_paired_R1=$( echo $TRIM_FILES_PAIRED_R1 | tr ":" "\n" | head -$sample_number | tail -1)
trimmedFastqArray_paired_R2=$( echo $TRIM_FILES_PAIRED_R2 | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                                                                                                                                                                                                                                                        
trimmedFastqArray_unpaired_R1=$( echo $TRIM_FILES_UNPAIRED_R1 | tr ":" "\n" | head -$sample_number | tail -1)
trimmedFastqArray_unpaired_R2=$( echo $TRIM_FILES_UNPAIRED_R2 | tr ":" "\n" | head -$sample_number | tail -1)   
TRIM_ARGS=$(echo $TRIM_ARGS | tr "_" " ")

echo -e "Running Trimmomatic for $SAMPLE....\n"                                                                                                                                                                                         
                                                                                                                                                                                                                                    
# Results folder per sample creation                                                                                                                                                                                               
mkdir -p $OUTPUT_DIR/QC/trimmomatic/$SAMPLE                                                                                                                                                                                             

java -jar $TRIMMOMATIC_PATH/trimmomatic-$trimmomatic_version.jar PE -phred33 $INPUT_DIR/$FASTQ_R1 $INPUT_DIR/$FASTQ_R2 $OUTPUT_DIR/QC/trimmomatic/$SAMPLE/$trimmedFastqArray_paired_R1 $OUTPUT_DIR/QC/trimmomatic/$SAMPLE/$trimmedFastqArray_unpaired_R1 $OUTPUT_DIR/QC/trimmomatic/$SAMPLE/$trimmedFastqArray_paired_R2 $OUTPUT_DIR/QC/trimmomatic/$SAMPLE/$trimmedFastqArray_unpaired_R2 $TRIM_ARGS
	                                                                                                                                                                                                                                                                                                                                                           
echo -e "Trimmomatic for $SAMPLE finished \n\n"                                                                                                                                                                                  
