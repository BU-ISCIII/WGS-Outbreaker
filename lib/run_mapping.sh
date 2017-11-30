#!/bin/bash
## Author S. Monzon & A.Hernández
## version v2.0

if [ $# -eq 0  ]; then
	echo -e "\n1º Execute mapping with bwa mem"
	echo -e "2º Convert SAM to BAM with samtools"
	echo -e "3º Mark duplicates with picadtools\n"
        echo -e "Usage: run_mapping.sh <config.file>"
        exit
fi


#Execute processing_config.sh

CONFIG_FILE=$1

# Check if run_outbreak_wgs.sh was execute
if [ -z $SCRIPTS_DIR ]; then
        SCRIPTS_DIR=$( cat $CONFIG_FILE | grep -w 'SCRIPTS_DIR' | cut -d '=' -f2 )
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

# Or other runner was execute
else
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"
fi


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


## Folder creation
echo -e "Creating $output_dir/Alignment"
mkdir -p $output_dir/Alignment
mkdir -p $output_dir/Alignment/SAM
mkdir -p $output_dir/Alignment/BAM

jobid_compress=$(cat $output_dir/logs/jobids.txt | grep -w "COMPRESS_FILE" | cut -d ':' -f2 )

#Check if trimmomatic was executed for input files
if [ -d $output_dir/QC/trimmomatic ]; then
	input_dir=$output_dir/QC/trimmomatic
	if [ "$use_sge" = "1" ]; then
		mapping_args="${SGE_ARGS} -pe openmp $threads -hold_jid ${jobid_compress}"
	fi

	#Mapping command with trimming files
        mapping_cmd="$SCRIPTS_DIR/bwa.sh \
                $threads \
                $input_dir \
                $output_dir/Alignment/SAM \
                $samples \
                $compress_paired_R1_list \
                $compress_paired_R2_list \
                $mappingArray_sam_list \
                $ref_path"

else
  	input_dir=$output_dir/raw
	#Mapping command with raw files
        mapping_cmd="$SCRIPTS_DIR/bwa.sh \
                $threads \
                $input_dir \
                $output_dir/Alignment/SAM \
                $samples \
                $fastq_R1_list \
                $fastq_R2_list  \
                $mappingArray_sam_list \
                $ref_path"
fi

samtobam_cmd="$SCRIPTS_DIR/samTobam.sh \
	$output_dir/Alignment/SAM \
	$output_dir/Alignment/BAM \
	$samples \
	$mappingArray_sam_list \
	$mappingArray_bam_list \
	$mappingArray_sorted_list \
	$mappingArray_rg_list \
	$picard_path \
	$platform \
	$model \
	$date_run \
	$library \
	$sequencing_center \
	$run_platform"

#Execute mapping and samTObam
if [ $mapping == "YES" ]; then

 	#In HPC
	if [ "$use_sge" = "1" ]; then	
		mapping_qsub=$( qsub $mapping_args -t 1-$sample_count -N $JOBNAME.MAPPING_BWA  $mapping_cmd)
		jobid_mapping=$( echo $mapping_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "MAPPING_BWAmem:$jobid_mapping\n" >> $output_dir/logs/jobids.txt
		
		samtobam_args="${SGE_ARGS} -pe openmp $threads -l h_vmem=$vmem -hold_jid $jobid_mapping"
		samtobam=$( qsub $samtobam_args -t 1-$sample_count -N $JOBNAME.SAMTOBAM $samtobam_cmd)
		jobid_samtobam=$( echo $samtobam | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "SAMTOBAM:$jobid_samtobam\n" >> $output_dir/logs/jobids.txt
	
	#In local
	else                              
		for count in `seq 1 $sample_count`
		do
			echo "Running mapping on sample $count"
			execute_mapping=$($mapping_cmd $count)
      			samtobam=$($samtobam_cmd $count)
     		done
	fi
fi

picard_cmd="$SCRIPTS_DIR/picard_duplicates.sh \
	$output_dir/Alignment/BAM \
	$samples \
	$mappingArray_rg_list \
	$duplicateBamArray_list \
	$picard_path"

#Execute duplicate filter
if [ $duplicate_filter == "YES" ]; then
	#In HPC
	if [ "$use_sge" = "1" ]; then
		picard_args="${SGE_ARGS} -pe openmp $threads -l h_vmem=$vmem -hold_jid $jobid_samtobam"
   		picard=$( qsub $picard_args -t 1-$sample_count -N $JOBNAME.PICARD_DUPLICATES $picard_cmd)
		jobid_picard=$( echo $picard | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "PICARD_DUPLICATES:$jobid_picard\n" >> $output_dir/logs/jobids.txt

	#Or local
  	else
		for count in `seq 1 $sample_count`
		do
			echo "Running mapping on sample $count"
			picard=$( $picard_cmd $count)
		done
	fi
fi
