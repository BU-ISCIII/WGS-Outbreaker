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
                                                                                                                                                                                                                                      
if [ $# != 14 -a "$use_sge" == "1" ]; then                                                                                                                                                                                           
  	echo "usage: ............"                                                                                                                                                                                                       
  	exit                                                                                                                                                                                                                             
elif [ $# != 15 -a "$use_sge" == "0" ]; then                                                                                                                                                                                         
  	echo "usage: ............"                                                                                                                                                                                                       
   	exit                                                                                                                                                                                                                             
fi                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                      
#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed                      
set -x                                                                                                                                                                                                                              
echo `date`                                                                                                                                                                                                                         
                                                                                                                                                                                                                                      
# Variables                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
DIR_SAM=$1
DIR_BAM=$2
SAMPLE_NAMES=$3
OUTPUT_SAM_NAMES=$4
OUTPUT_BAM_NAMES=$5
OUTPUT_BAM_SORTED_NAMES=$6
OUTPUT_BAM_SORTED_RG_NAMES=$7
PLATFORM=$8
MODEL=$9
DATE_RUN=${10}
LIBRARY=${11}
SEQUENCING_CENTER=${12}
RUN_PLATFORM=${13}
PICARD_PATH=${14}

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                                        
  	sample_number=$SGE_TASK_ID                                                                                                                                                                                                       
else                                                                                                                                                                                                                                 
  	sample_number=${15}                                                                                                                                                                                                                 
fi                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                      
SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                        
OUTPUT_SAM_NAME=$( echo $OUTPUT_SAM_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                         
OUTPUT_BAM_NAME=$( echo $OUTPUT_BAM_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                          
OUTPUT_BAM_SORTED_NAME=$( echo $OUTPUT_BAM_SORTED_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                          
OUTPUT_BAM_SORTED_RG_NAME=$( echo $OUTPUT_BAM_SORTED_RG_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                           

mkdir -p $DIR_BAM/$SAMPLE

samtools view -bS $DIR_SAM/$SAMPLE/$OUTPUT_SAM_NAME -o $DIR_BAM/$SAMPLE/$OUTPUT_BAM_NAME                         
samtools sort -o $DIR_BAM/$SAMPLE/$OUTPUT_BAM_SORTED_NAME -T $DIR_BAM/$SAMPLE/$OUTPUT_BAM_SORTED_NAME $DIR_BAM/$SAMPLE/$OUTPUT_BAM_NAME                           
                                                                                                                                                     
java $JAVA_RAM -jar $PICARD_PATH/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT INPUT=$DIR_BAM/$SAMPLE/$OUTPUT_BAM_SORTED_NAME OUTPUT=$DIR_BAM/$SAMPLE/$OUTPUT_BAM_SORTED_RG_NAME RGID=$DATE_RUN-$LIBRARY-$MODEL-$PLATFORM-$SEQUENCING_CENTER RGLB=$LIBRARY RGPL=$PLATFORM RGSM=$SAMPLE RGPU=$RUN_PLATFORM RGDT=$DATE_RUN RGCN=$SEQUENCING_CENTER

rm -r $DIR_SAM/$SAMPLE                                                                                                           
                                                                                                                                                     
# Se indexa el fichero BAM con samtools.                                                                                                           
samtools index $DIR_BAM/$SAMPLE/$OUTPUT_BAM_SORTED_RG_NAME                                                                
