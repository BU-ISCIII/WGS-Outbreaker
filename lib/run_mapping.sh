#!/bin/bash
## Author S. Monzon
## version v2.0

# Help 
# usage: align.sh ....
#

set -e
set -u
set -x

# Variables
MAPPING=$1
USE_SGE=$2
OUTPUT_DIR=$3
REF_PATH=$4
THREADS=$5
FASTQ_R1_list=$6
FASTQ_R2_list=$7
sample_number=$8
SAMPLE_NAMES=$9
OUTPUT_SAM_NAMES=${10}
OUTPUT_BAM_NAMES=${11}
OUTPUT_BAM_SORTED_NAMES=${12}
BAMSTAT_PRE_NAMES=${13}
BAMSTAT_POST_NAMES=${14}
OUTPUT_DUPLICATE_NAMES=${15}
DUPLICATE_FILTER=${16}
PICARD_PATH=${17}
EXOME_ENRICHMENT=${18}
PLATFORM=${19}
MODEL=${20}
DATE_RUN=${21}
LIBRARY=${22}
SEQUENCING_CENTER=${23}
RUN_PLATFORM=${24}
OUTPUT_BAM_SORTED_RG_NAMES=${25}

## Folder creation
echo -e "Creating $OUTPUT_DIR/Alignment"
mkdir -p $OUTPUT_DIR/Alignment
mkdir -p $OUTPUT_DIR/Alignment/SAM
mkdir -p $OUTPUT_DIR/Alignment/BAM

jobid_trimmomatic=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )

if [ -d $OUTPUT_DIR/QC/trimmomatic ]; then
	DIR=$OUTPUT_DIR/QC/trimmomatic
	MAPPING_ARGS="${SGE_ARGS} -pe orte $THREADS -hold_jid ${jobid_trimmomatic}"
else
  	DIR=$OUTPUT_DIR/raw
fi

MAPPING_CMD="$SCRIPTS_DIR/bwa.sh $DIR $THREADS $REF_PATH $OUTPUT_DIR/Alignment/SAM $FASTQ_R1_list $FASTQ_R2_list $SAMPLE_NAMES $OUTPUT_SAM_NAMES"  
SAMTOBAM_CMD="$SCRIPTS_DIR/samTobam.sh $OUTPUT_DIR/Alignment/SAM $OUTPUT_DIR/Alignment/BAM $SAMPLE_NAMES $OUTPUT_SAM_NAMES $OUTPUT_BAM_NAMES $OUTPUT_BAM_SORTED_NAMES $OUTPUT_BAM_SORTED_RG_NAMES $PLATFORM $MODEL $DATE_RUN $LIBRARY $SEQUENCING_CENTER $RUN_PLATFORM $PICARD_PATH"

if [ $MAPPING == "YES" ]; then                                                                                                                                                                                                                                                                                       
  	if [ "$USE_SGE" = "1" ]; then
  		MAPPING=$( qsub $MAPPING_ARGS -t 1-$sample_number -N $JOBNAME.MAPPING $MAPPING_CMD)
     	jobid_mapping=$( echo $MAPPING | cut -d ' ' -f3 | cut -d '.' -f1 )                                                                                                                                                                                                                                    
     	echo -e "MAPPING:$jobid_mapping\n" >> $OUTPUT_DIR/logs/jobids.txt

  		SAMTOBAM_ARGS="$SGE_ARGS -hold_jid $jobid_mapping"
  		SAMTOBAM=$( qsub $SAMTOBAM_ARGS -t 1-$sample_number -N $JOBNAME.SAMTOBAM $SAMTOBAM_CMD) 
      	jobid_samtobam=$( echo $SAMTOBAM | cut -d ' ' -f3 | cut -d '.' -f1 )                                                                                                                                                                                                                                     
      	echo -e "SAMTOBAM:$jobid_samtobam\n" >> $OUTPUT_DIR/logs/jobids.txt

 	else                                                                                                                                                                                                                                                                                                              
      	for count in `seq 1 $sample_number`                                                                                                                                                                                                                                                                           
      	do                                                                                                                                                                                                                                                                                                            
      		echo "Running mapping on sample $count"                                                                                                                                                                                                                                                               
      		MAPPING=$($MAPPING_CMD $count)
      		SAMTOBAM=$($SAMTOBAM_CMD $count)
     	done                                                                                                                                                                                                                                                                                                          
    fi                                                                                                                                                                                                                                                                                                                
fi

PICARD_CMD="$SCRIPTS_DIR/picard_duplicates.sh $OUTPUT_DIR/Alignment/BAM $SAMPLE_NAMES $OUTPUT_BAM_SORTED_RG_NAMES $OUTPUT_DUPLICATE_NAMES $PICARD_PATH"
if [ $DUPLICATE_FILTER == "YES" ]; then                                                                                                                                                                                                                                                                                                         
   	if [ "$USE_SGE" = "1" ]; then                                                                                                                                                                                                                                                                                                      
   		PICARD_ARGS="$SGE_ARGS -hold_jid $jobid_samtobam"
   		PICARD=$( qsub $PICARD_ARGS -t 1-$sample_number -N $JOBNAME.PICARD $PICARD_CMD)                                                                                                                                                                                                                                                
      	jobid_picard=$( echo $PICARD | cut -d ' ' -f3 | cut -d '.' -f1 )                                                                                                                                                                                                                                                             
      	echo -e "PICARD:$jobid_picard\n" >> $OUTPUT_DIR/logs/jobids.txt
  	else                                                                                                                                                                                                                                                                                                                               
       	for count in `seq 1 $sample_number`                                                                                                                                                                                                                                                                                            
       	do                                                                                                                                                                                                                                                                                                                             
       		echo "Running mapping on sample $count"                                                                                                                                                                                                                                                                                    
       		PICARD=$($PICARD_CMD $count)                                                                                                                                                                                                                                                                                             
      	done                                                                                                                                                                                                                                                                                                                           
     fi                                                                                                                                                                                                                                                                                                                                 
 fi

 ## FASTQC                                                                                                                                                                      
 BAMUTIL_PREDUPLICATES="$SCRIPTS_DIR/bamutil.sh $OUTPUT_DIR/Alignment/BAM $SAMPLE_NAMES $OUTPUT_BAM_SORTED_RG_NAMES $BAMSTAT_PRE_NAMES $EXOME_ENRICHMENT"                                                   
 BAMUTIL_POSTDUPLICATES="$SCRIPTS_DIR/bamutil.sh $OUTPUT_DIR/Alignment/BAM $SAMPLE_NAMES $OUTPUT_DUPLICATE_NAMES $BAMSTAT_POST_NAMES $EXOME_ENRICHMENT"  
                                                                                                                                                                                
 if [ "$USE_SGE" = "1" ]; then                                                                                                                                                  
 	    BAMUTIL_PRE_ARGS="${SGE_ARGS} -hold_jid $jobid_samtobam"
 	 	BAMUTIL_PRE=$( qsub $BAMUTIL_PRE_ARGS -t 1-$sample_number -N $JOBNAME.BAMUTIL_PRE $BAMUTIL_PREDUPLICATES)                                                                                
        jobid_bamutil_pre=$( echo $BAMUTIL_PRE | cut -d ' ' -f3 | cut -d '.' -f1 )                                                                                                   
        echo -e "BAMUTIL_PRE:$jobid_bamutil_pre\n" >> $OUTPUT_DIR/logs/jobids.txt
 	if [ $DUPLICATE_FILTER == "YES" ]; then                                                                                                                                            
 		BAMUTIL_POST_ARGS="${SGE_ARGS} -hold_jid $jobid_picard"                                                                                                                 
 		BAMUTIL_POST=$( qsub $BAMUTIL_POST_ARGS -t 1-$sample_number -N $JOBNAME.BAMUTIL_POST $BAMUTIL_POSTDUPLICATES)                                                                       
 		jobid_bamutil_post=$( echo $BAMUTIL_POST | cut -d ' ' -f3 | cut -d '.' -f1 )
 		echo -e "BAMUTIL_POST:$jobid_bamutil_post\n" >> $OUTPUT_DIR/logs/jobids.txt
 	fi                                                                                                                                                                         
 else                                                                                                                                                                           
     for count in `seq 1 $sample_number`                                                                                                                                        
     do                                                                                                                                                                         
     	echo "Running fastqc on sample $count"                                                                                                                                 
     	BAMUTIL_PRE=$($BAMUTIL_PREDUPLICATES $count)                                                                                                                               
 		if [ $DUPLICATE_FILTER == "YES" ]; then                                                                                                                                        
 			BAMUTIL_POST=$( $BAMUTIL_POSTDUPLICATES $count)                                                                                                                         
 		fi                                                                                                                                                                     
     done                                                                                                                                                                       
 fi                                                                                                                                                                             
