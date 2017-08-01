 #!/bin/bash                                          
 ## Author S. Monzon                                  
 ## version v2.0                                      
                                                      
 # Help
 # usage: run_preprocessing.sh <INPUT_DIR> <RESULTS_DIR> <SAMPLES> 
 # 

 ## VARIABLES
set -e
set -u
set -x

TRIMMING=$1
USE_SGE=$2 
INPUT_DIR=$3 
OUTPUT_DIR=$4
THREADS=$5
fastq_R1_list=$6 
fastq_R2_list=$7
sample_number=$8
sample_names=$9
TRIM_ARGS=${10}
trimmomatic_version=${11}
TRIMMOMATIC_PATH=${12}
trimmedFastqArray_paired_R1_list=${13}
trimmedFastqArray_paired_R2_list=${14}
trimmedFastqArray_unpaired_R1_list=${15}
trimmedFastqArray_unpaired_R2_list=${16}

## Create directories
date                            
echo -e "Creating $OUTPUT_DIR" 
mkdir -p $OUTPUT_DIR/QC          
mkdir -p $OUTPUT_DIR/QC/fastqc    


## TRIMMOMATIC

TRIMMING_CMD="$SCRIPTS_DIR/trimmomatic.sh $INPUT_DIR $OUTPUT_DIR $sample_names $fastq_R1_list $fastq_R2_list $trimmedFastqArray_paired_R1_list $trimmedFastqArray_paired_R2_list $trimmedFastqArray_unpaired_R1_list $trimmedFastqArray_unpaired_R2_list $THREADS $TRIM_ARGS $trimmomatic_version $TRIMMOMATIC_PATH" 

if [ $TRIMMING == "YES" ]; then
	mkdir -p $OUTPUT_DIR/QC/trimmomatic
 	if [ "$USE_SGE" = "1" ]; then                                                                                           
 		TRIMMOMATIC=$( qsub $SGE_ARGS -pe orte $THREADS -t 1-$sample_number -N $JOBNAME.TRIMMOMATIC $TRIMMING_CMD)                                 
    	jobid_trimmomatic=$( echo $TRIMMOMATIC | cut -d ' ' -f3 | cut -d '.' -f1 )                                                    
    	echo -e "TRIMMOMATIC:$jobid_trimmomatic\n" >> $OUTPUT_DIR/logs/jobids.txt
	else                                                                                                                    
     	for count in `seq 1 $sample_number`                                                                                 
     	do                                                                                                                  
     		echo "Running trimmomatic on sample $count"                                                                          
     		TRIMMOMATIC=$($TRIMMING_CMD $count)                                                                            
    	done                                                                                                                
    fi                                                                                                                      
fi

## FASTQC
FASTQC_PRETRIMMING="$SCRIPTS_DIR/fastqc.sh $OUTPUT_DIR/raw $OUTPUT_DIR $sample_names $fastq_R1_list $fastq_R2_list $THREADS"
FASTQC_POSTRIMMING="$SCRIPTS_DIR/fastqc.sh $OUTPUT_DIR/QC/trimmomatic $OUTPUT_DIR $sample_names $trimmedFastqArray_paired_R1_list $trimmedFastqArray_paired_R2_list $THREADS" 

if [ "$USE_SGE" = "1" ]; then
	FASTQC_PRE=$( qsub $SGE_ARGS -pe orte $THREADS -t 1-$sample_number -N $JOBNAME.FASTQ_PRE $FASTQC_PRETRIMMING)
    jobid_fastqc_pre=$( echo $FASTQC_PRE | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "FASTQC_PRE:$jobid_fastqc_pre\n" >> $OUTPUT_DIR/logs/jobids.txt 
	if [ $TRIMMING == "YES" ]; then
		FASTQC_ARGS="${SGE_ARGS} -pe orte $THREADS -hold_jid $jobid_trimmomatic"
		FASTQC_POST=$( qsub $FASTQC_ARGS -t 1-$sample_number -N $JOBNAME.FASTQ_POST $FASTQC_POSTRIMMING)
		jobid_fastqc_post=$( echo $FASTQC_POST | cut -d ' ' -f3 | cut -d '.' -f1 ) 
 		echo -e "FASTQC_POST:$jobid_fastqc_post\n" >> $OUTPUT_DIR/logs/jobids.txt  
	fi
else
    for count in `seq 1 $sample_number`
    do 
    	echo "Running fastqc on sample $count"
    	FASTQC_PRE=$($FASTQC_PRETRIMMING $count)
		if [ $TRIMMING == "YES" ]; then
			FASTQC_POST=$( $FASTQC_POSTRIMMING $count)
		fi
    done
fi
