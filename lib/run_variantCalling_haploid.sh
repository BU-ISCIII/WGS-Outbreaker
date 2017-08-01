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
DUPLICATES=$3
OUTPUT_DIR=$4
REF_PATH=$5
THREADS=$6
SAMPLE_NAMES=$7
OUTPUT_BAM_NAMES=$8
VCF_NAMES=$9
sample_number=${10}
OUTPUT_BAM_REALIGNED_NAMES=${11}
OUTPUT_BAM_RECALIBRATED_NAMES=${12}
GATK_PATH=${13}
SNPS_NAME=${14}
SNPS_FIL_NAME=${15}
INDELS_NAME=${16}
INDELS_FIL_NAME=${17}
VCF_FIL_NAME=${18}
KNOWN_SNPS=${19}
KNOWN_INDELS=${20}

## Folder creation
echo -e "Creating $OUTPUT_DIR/variant_calling"
mkdir -p $OUTPUT_DIR/variant_calling

if [ "$USE_SGE" = "1" -a $DUPLICATES == "YES" ]; then
 	jobid=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "PICARD" | cut -d ':' -f2 )
 	PRECALLING_ARGS="${SGE_ARGS} -pe orte $THREADS -hold_jid $jobid"
else
 	jobid=$(cat $OUTPUT_DIR/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )
 	PRECALLING_ARGS="${SGE_ARGS} -pe orte $THREADS -hold_jid $jobid"
                                                                                                                                                                                                                                                                                                                                                                                                  
fi

mkdir -p $OUTPUT_DIR/variant_calling/variants_gatk
PRECALLING_CMD="$SCRIPTS_DIR/gatk_preprocessing.sh $OUTPUT_DIR/Alignment/BAM $THREADS $REF_PATH $SAMPLE_NAMES $KNOWN_SNPS $KNOWN_INDELS $OUTPUT_BAM_NAMES $OUTPUT_DIR/variant_calling/variants_gatk $OUTPUT_BAM_REALIGNED_NAMES $OUTPUT_BAM_RECALIBRATED_NAMES $GATK_PATH"

if [ $KNOWN_SNPS == "NO" ];then
	CALLING_CMD="$SCRIPTS_DIR/gatk_haploid.sh $OUTPUT_DIR/variant_calling/variants_gatk/realignment $THREADS $REF_PATH $OUTPUT_DIR/variant_calling/variants_gatk $SAMPLE_NAMES $OUTPUT_BAM_REALIGNED_NAMES $VCF_NAMES $SNPS_NAME $SNPS_FIL_NAME $INDELS_NAME $INDELS_FIL_NAME $VCF_FIL_NAME $GATK_PATH"
else
  	CALLING_CMD="$SCRIPTS_DIR/gatk_haploid.sh $OUTPUT_DIR/variant_calling/variants_gatk/recalibration $THREADS $REF_PATH $OUTPUT_DIR/variant_calling/variants_gatk $SAMPLE_NAMES $OUTPUT_BAM_RECALIBRATED_NAMES $VCF_NAMES $SNPS_NAME $SNPS_FIL_NAME $INDELS_NAME $INDELS_FIL_NAME $VCF_FIL_NAME $GATK_PATH" 
fi

if [ $VARIANT_CALLING == "YES" ]; then
		if [ "$USE_SGE" = "1" ]; then
            PRECALLING=$( qsub $PRECALLING_ARGS -t 1-$sample_number -N $JOBNAME.CALLING $PRECALLING_CMD)
    		jobid_precalling=$( echo $PRECALLING | cut -d ' ' -f3 | cut -d '.' -f1 )
    		CALLING_ARGS="${SGE_ARGS} -hold_jid $jobid_precalling"
    		CALLING=$( qsub $CALLING_ARGS -N $JOBNAME.CALLING $CALLING_CMD)
       		jobid_calling=$( echo $CALLING | cut -d ' ' -f3 | cut -d '.' -f1 )
       		echo -e "Variant Calling:$jobid_precalling - $jobid_calling \n" >> $OUTPUT_DIR/logs/jobids.txt
		else
        	for count in `seq 1 $sample_number`
        	do
        		echo "Running variant calling on sample $count"
        		PRECALLING=$($PRECALLING_CMD $count)
       		done
        	CALLING=$($CALLING_CMD)
      fi
fi
