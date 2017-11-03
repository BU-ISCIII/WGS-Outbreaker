 #!/bin/bash                                          
 ## Author S. Monzon                                  
 ## version v2.0                                      
                                                      
 # Help
 # usage: run_preprocessing.sh <INPUT_DIR> <RESULTS_DIR> <SAMPLES> 

CONFIG_FILE=$1

#Execute processing_config.sh
if [ -z $SCRIPTS_DIR ]; then
        SCRIPTS_DIR=$( cat $CONFIG_FILE | grep -w 'SCRIPTS_DIR' | cut -d '=' -f2 )
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

else
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"
fi


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


## Create directories
date                            
echo -e "Creating $output_dir" 
mkdir -p $output_dir/QC          
mkdir -p $output_dir/QC/fastqc    


## TRIMMOMATIC

trimming_cmd="$SCRIPTS_DIR/trimmomatic.sh \
	$threads\
	$input_dir \
	$output_dir \
	$samples \
	$fastq_R1_list \
	$fastq_R2_list \
	$trimmedFastqArray_paired_R1_list \
	$trimmedFastqArray_paired_R2_list \
	$trimmedFastqArray_unpaired_R1_list \
	$trimmedFastqArray_unpaired_R2_list \
	$trim_args \
	$trimmomatic_version \
	$trimmomatic_path" 

if [ $trimming == "YES" ]; then
	mkdir -p $output_dir/QC/trimmomatic
 	if [ "$use_sge" = "1" ]; then                                                                                           
 		trimmomatic_arg="${SGE_ARGS} -pe openmp $threads -l h_vmem=$vmem -t 1-$sample_count"
		trimmomatic=$( qsub $trimmomatic_arg -N $JOBNAME.TRIMMOMATIC $trimming_cmd) 
		jobid_trimmomatic=$( echo $trimmomatic | cut -d ' ' -f3 | cut -d '.' -f1 )                                              
    		echo -e "TRIMMOMATIC:$jobid_trimmomatic\n" >> $output_dir/logs/jobids.txt
	else                                                                                                                    
     		for count in `seq 1 $sample_count`
		do
			echo "Running trimmomatic on sample $count"
			trimmomatic=$($trimming_cmd $count)
		done                                                                                                                
	fi
fi

## FASTQC
fastqc_pretrimming_cmd="$SCRIPTS_DIR/fastqc.sh \
	$threads \
	$output_dir/raw \
	$output_dir \
	$samples \
	$fastq_R1_list \
	$fastq_R2_list"

fastqc_postrimming_cmd="$SCRIPTS_DIR/fastqc.sh \
	$threads \
	$output_dir/QC/trimmomatic \
	$output_dir \
	$samples \
	$trimmedFastqArray_paired_R1_list \
	$trimmedFastqArray_paired_R2_list " 

##COMPRESSING FILES
compress_file_cmd="$SCRIPTS_DIR/compress.sh \
        $output_dir/QC/trimmomatic \
        $output_dir/CQ/trimmomatic \
        $samples \
        $trimmedFastqArray_paired_R1_list \
        $trimmedFastqArray_paired_R2_list \
        $trimmedFastqArray_unpaired_R1_list \
        $trimmedFastqArray_unpaired_R2_list"


if [ "$use_sge" = "1" ]; then
	fastqc_pre=$( qsub $SGE_ARGS -pe openmp $threads -t 1-$sample_count -N $JOBNAME.FASTQ_PRE $fastqc_pretrimming_cmd)
    jobid_fastqc_pre=$( echo $fastqc_pre | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "FASTQC_PRE:$jobid_fastqc_pre\n" >> $output_dir/logs/jobids.txt 
	if [ $trimming == "YES" ]; then
		fastqc_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_trimmomatic"
		fastqc_post=$( qsub $fastqc_arg -t 1-$sample_count -N $JOBNAME.FASTQ_POST $fastqc_postrimming_cmd)
		jobid_fastqc_post=$( echo $fastqc_post | cut -d ' ' -f3 | cut -d '.' -f1 ) 
 		echo -e "FASTQC_POST:$jobid_fastqc_post\n" >> $output_dir/logs/jobids.txt 

		compress_arg="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_fastqc_post"
                compress_file=$( qsub $compress_arg -t 1-$sample_count -N $JOBNAME.COMPRESS_FILE $compress_file_cmd)
                jobid_compress=$( echo $compress_file | cut -d ' ' -f3 | cut -d '.' -f1 ) 
		echo -e "COMPRESS_FILE:$jobid_compress\n" >> $output_dir/logs/jobids.txt
	fi
else
    for count in `seq 1 $sample_count`
    do 
    	echo "Running fastqc on sample $count"
    	fastqc_pre=$($fastqc_pretrimming_cmd $count)
		if [ $trimming == "YES" ]; then
			fastqc_post=$($fastqc_postrimming_cmd $count)
			compress_file=$($compress_file_cmd $count)
		fi
    done
fi
