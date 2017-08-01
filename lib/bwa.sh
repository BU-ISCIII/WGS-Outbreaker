#!/bin/bash                                                                                                                                                                                                                        
## Author S. Monzon                                                                                                                                                                                                                
## version v2.0                                                                                                                                                                                                                    
                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                    
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
                                                                                                                                                                                                                                    
if [ $# != 8 -a "$use_sge" == "1" ]; then                                                                                                                                                                                          
 	echo "usage: ............"                                                                                      
 	exit                                                                                                                                                                                                                           
elif [ $# != 9 -a "$use_sge" == "0" ]; then                                                                                                                                                                                        
 	echo "usage: ............"                                                                             
  	exit                                                                                                                                                                                                                           
fi                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                    
#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed                     
set -x                                                                                                                                                                                                                             
echo `date` 

# Variables

DIR=$1
THREADS=$2
REF_PATH=$3
OUTPUT_DIR=$4
FASTQ_FILES_R1=$5
FASTQ_FILES_R2=$6
SAMPLE_NAMES=$7
OUTPUT_SAM_NAMES=$8

if [ "$use_sge" = "1" ]; then                                                      
 	sample_number=$SGE_TASK_ID                                                     
else                                                                               
 	sample_number=${9}                                                               
fi                                                                                 
                                                                                                                                                                        
SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)       
FASTQ_R1=$( echo $FASTQ_FILES_R1 | tr ":" "\n" | head -$sample_number | tail -1)   
FASTQ_R2=$( echo $FASTQ_FILES_R2 | tr ":" "\n" | head -$sample_number | tail -1)   
OUTPUT_SAM_NAME=$( echo $OUTPUT_SAM_NAMES | tr ":" "\n" | head -$sample_number | tail -1)    

mkdir -p $OUTPUT_DIR/$SAMPLE

echo -e "2) Alignment using BWA: $SAMPLE"                                                                                                                                                            
                                                                                                                                                                                                       
echo "Alignment PE"                                                                                                                                                                            
bwa aln -t $THREADS $REF_PATH $DIR/$SAMPLE/$FASTQ_R1 > $TEMP/$SAMPLE-align1.sai                                                                                                                                        
bwa aln -t $THREADS $REF_PATH $DIR/$SAMPLE/$FASTQ_R2 > $TEMP/$SAMPLE-align2.sai                                                                                                                                        
                                                                                                                                                                                                       
 # Use sampe with paired end data.                                                                                                                                                                    
bwa sampe -P $REF_PATH $TEMP/$SAMPLE-align1.sai $TEMP/$SAMPLE-align2.sai $DIR/$SAMPLE/$FASTQ_R1 $DIR/$SAMPLE/$FASTQ_R2 > $OUTPUT_DIR/$SAMPLE/$OUTPUT_SAM_NAME   
