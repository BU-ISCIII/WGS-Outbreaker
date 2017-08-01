#!/bin/bash                                                                                                                                                                                                                                                                                                            
## Author S. Monzon                                                                                                                                                                                                                                                                                                    
## version v2.0                                                                                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                                                                                                        
# Help                                                                                                                                                                                                                                                                                                                 
# usage: run_variantCalling.sh ....                                                                                                                                                                                                                                                                                                 
#                                                                                                                                                                                                                                                                                                                      
                                                                                                                                                                                                                                                                                                                        
set -e                                                                                                                                                                                                                                                                                                                 
set -u                                                                                                                                                                                                                                                                                                                 
set -x                                                                                                                                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                                                        
# Variables
USE_SGE=$1
VARIANT_CALLING=$2
VARIANT_CALLER=$3
DUPLICATES=$4
OUTPUT_DIR=$5                                                                                                                                                                                                                                                                                                         
REF_PATH=$6                                                                                                                                                                                                                                                                                                            
THREADS=$7                                                                                                                                                                                                                                                                                                             
SAMPLE_NAMES=$8                                                                                                                                                                                                                                                                                                        
OUTPUT_BAM_NAMES=$9                                                                                                                                                                                                                                                                                                 
EXOME_ENRICHMENT=${10}
CONTROL_NAMES=${11}
CASE_NAMES=${12}
VCF_NAMES=${13}
STRELKA_CONFIG=${14}
sample_number=${15}
                                                                                                                                                                                                                                                                                                                        
## Folder creation                                                                                                                                                                                                                                                                                                     
echo -e "Creating $OUTPUT_DIR/variant_calling"                                                                                                                                                                                                                                                                               
mkdir -p $OUTPUT_DIR/variant_calling

if [ $DUPLICATES == "YES" ]; then
	jobid=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "PICARD" | cut -d ':' -f2 )
	CALLING_ARGS="${SGE_ARGS} -pe orte $THREADS -hold_jid $jobid"
else
	jobid=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )
	CALLING_ARGS="${SGE_ARGS} -pe orte $THREADS -hold_jid $jobid" 
fi

if [ $VARIAN_CALLING == "YES" ];then
	if [ $VARIANT_CALLER == "STRELKA" ];then
		mkdir -p $OUTPUT_DIR/variant_calling/variants_strelka                                                                                                                                                                                                                                                                                                    
		CALLING_CMD="$SCRIPTS_DIR/strelka.sh $OUTPUT_DIR/Alignment/BAM $THREADS $REF_PATH $OUTPUT_DIR/variant_calling/variants_strelka $SAMPLE_NAMES $OUTPUT_BAM_NAMES $CONTROL_NAMES $CASE_NAMES $VCF_NAMES $STRELKA_CONFIG"                                                                                                                                                                       
		if [ "$USE_SGE" = 1 ];then
           calling_sge $CALLING_CMD $CALLING_ARGS $sample_number
		else
           calling $CALLING_CMD $sample_number
		fi
	fi
fi


if [ $VARIANT_CALLING == "YES" ]; then                                                                                                                                                                                                                                                                                         


fi

function calling_sge  {
    	CALLING=$( qsub $2 -t 1-$3 -N $JOBNAME.CALLING $1)                                                                                                                                                                                                                                     
       	jobid_calling=$( echo $CALLING | cut -d ' ' -f3 | cut -d '.' -f1 )                                                                                                                                                                                                                                                      
       	echo -e "Variant Calling:$jobid_calling\n" >> $OUTPUT_DIR/logs/jobids.txt                                                                                                                                                                                                                                                      	
}
function calling {
	for count in `seq 1 $2`                                                                                                                                                                                                                                                                                      
    do                                                                                                                                                                                                                                                                                                                       
    	CALLING=$($1 $count)                                                                                                                                                                                                                                                                                       
    done                                                                                                                                                                                                                                                                                                                     
}
