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
set -x

## Usage                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
if [ $# != 11 -a "$use_sge" == "1" ]; then                                                                                                                                                                                                                                                                                                                                                                                                                                               
   	echo "usage: ............"                                                                                                                                                                                                                                                                                                                                                                                                                                                          
   	exit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
elif [ $# != 12 -a "$use_sge" == "0" ]; then                                                                                                                                                                                                                                                                                                                                                                                                                                             
   	echo "usage: ............"                                                                                                                                                                                                                                                                                                                                                                                                                                                          
    	exit                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
fi                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed                                                                                                                                                                                                                                                                          
set -x                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
echo `date`                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
# Variables                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
DIR_BAM=$1 
THREADS=$2
REF_PATH=$3
SAMPLE_NAMES=$4
KNOWN_SNPS=$5
KNOWN_INDELS=$6
OUTPUT_BAM_NAMES=$7                                                                                                                                                                                                                                                                                                                                                                                                                                                              
OUTPUT_DIR=$8
OUTPUT_BAM_REALIGNED_NAMES=$9
OUTPUT_BAM_RECALIBRATED_NAMES=${10}
GATK_PATH=${11}

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                                                                                                                                                                                                                                                                                           
   	sample_number=$SGE_TASK_ID                                                                                                                                                                                                                                                                                                                                                                                                                                                          
else                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
   	sample_number=${12}                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
fi                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                                                                                                                                                                                                                                                                            
BAM_NAME=$( echo $OUTPUT_BAM_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                                                                                                                                                                                                                                                               
BAM_REALIGNED_NAME=$( echo $OUTPUT_BAM_REALIGNED_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                                                                                                                                                                                                                                                                
BAM_RECALIBRATED_NAME=$( echo $OUTPUT_BAM_RECALIBRATED_NAMES | tr ":" "\n" | head -$sample_number | tail -1)                                                                                                                                                                                                                                                                                                                                                                                                 

mkdir -p $OUTPUT_DIR/realignment

if [ $KNOWN_INDELS == "NO" ];then

	java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
 		-T RealignerTargetCreator \
		-I $DIR_BAM/$SAMPLE/$BAM_NAME \
 		-R $REF_PATH \
 		-o $OUTPUT_DIR/realignment/$BAM_NAME-IndelRealigner.intervals \
 		-nt $THREADS \
 		-S LENIENT \
 		-log $OUTPUT_DIR/realignment/$BAM_NAME-targetCreator.log
                                                                                       
	java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
 		-T IndelRealigner \
		-I $DIR_BAM/$SAMPLE/$BAM_NAME \
 		-R $REF_PATH \
 		-targetIntervals $OUTPUT_DIR/realignment/$BAM_NAME-IndelRealigner.intervals \
 		-o $OUTPUT_DIR/realignment/$BAM_REALIGNED_NAME \
 		-S LENIENT \
 		-log $OUTPUT_DIR/$BAM_NAME-realigner.log

else
    
    java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
    	-T RealignerTargetCreator \
    	-I $DIR_BAM/$SAMPLE/$BAM_NAME \
    	-known $KNOWN_INDELS \
    	-R $REF_PATH \
    	-o $OUTPUT_DIR/realignment/$BAM_NAME-IndelRealigner.intervals \
    	-nt $THREADS \
    	-S LENIENT \
    	-log $OUTPUT_DIR/realignment/$BAM_NAME-targetCreator.log
                                                                                        
    java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
    	-T IndelRealigner \
    	-I $DIR_BAM/$SAMPLE/$BAM_NAME \
    	-known $KNOWN_INDELS \
    	-R $REF_PATH \
    	-targetIntervals $OUTPUT_DIR/realignment/$BAM_NAME-IndelRealigner.intervals \
    	-o $OUTPUT_DIR/realignment/$BAM_REALIGNED_NAME \
    	-S LENIENT \
    	-log $OUTPUT_DIR/$BAM_NAME-realigner.log                                        

fi


if [ $KNOWN_SNPS != "NO" ];then
	echo -e "2) Base Quality Recalibration"                                               
                                                                                       
	mkdir -p $OUTPUT_DIR/recalibration
	
	java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
      	-T BaseRecalibrator \
      	-R $REF_PATH \
 	  	-I $OUTPUT_DIR/realignment/$BAM_REALIGNED_NAME \
 	  	-knownSites $KNOWN_SNPS \
 	  	-cov ReadGroupCovariate \
 	  	-cov QualityScoreCovariate \
 	  	-cov CycleCovariate \
 	  	-cov ContextCovariate \
      	-o $OUTPUT_DIR/recalibration/$BAM_NAME-recal1_data.grp \
 	  	-nct $THREADS \
 	  	-S LENIENT \
 	  	-log $OUTPUT_DIR/gatk_trio/$BAM_NAME-recal.log
                                                                                       
	java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
    	-T BaseRecalibrator \
    	-BQSR $OUTPUT_DIR/recalibration/$BAM_NAME-recal1_data.grp \
    	-R $REF_PATH \
    	-I $OUTPUT_DIR/realignment/$BAM_REALIGNED_NAME \
    	-knownSites $KNOWN_SNPS \
   		-cov ReadGroupCovariate \
    	-cov QualityScoreCovariate \
    	-cov CycleCovariate \
    	-cov ContextCovariate \
    	-o $OUTPUT_DIR/recalibration/$BAM_NAME-recal2_data.grp \
    	-nct $THREADS \
    	-S LENIENT \
    	-log $OUTPUT_DIR/$BAM_NAME-recal.log                                        
                                                                                       
	java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
    	-T AnalyzeCovariates \
    	-R $REF_PATH \
    	-before $OUTPUT_DIR/recalibration/$BAM_NAME-recal1_data.grp \
    	-after $OUTPUT_DIR/recalibration/$BAM_NAME-recal2_data.grp \
    	-csv $OUTPUT_DIR/recalibration/$BAM_NAME-BQSR.csv \
    	-plots $OUTPUT_DIR/recalibration/$BAM_NAME-BQSR.pdf \
    	-log $OUTPUT_DIR/$BAM_NAME-analyzecovariates.log                       
                                                                                       
	java -Djava.io.tmpdir=$TEMP $JAVA_RAM -jar $GATK_PATH/GenomeAnalysisTK.jar \
    	-T PrintReads \
    	-R $REF_PATH \
    	-I $OUTPUT_DIR/realignment/$BAM_REALIGNED_NAME \
    	-BQSR $OUTPUT_DIR/recalibration/$BAM_NAME-recal1_data.grp \
    	-o $OUTPUT_DIR/recalibration/$BAM_RECALIBRATED_NAME \
    	-nct $THREADS \
    	-S LENIENT \
    	-log $OUTPUT_DIR/$BAM_NAME-print_recal.log
fi
