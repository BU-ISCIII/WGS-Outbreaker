#!/bin/bash
## Author S. Monzon
## version v2.0

# Help 
# usage: align.sh ....
#

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


#Execure processing_config.sh
source $SCRIPTS_DIR/processing_config.sh

## Folder creation
echo -e "Creating $output_dir/Alignment"
mkdir -p $output_dir/Alignment
mkdir -p $output_dir/Alignment/SAM
mkdir -p $output_dir/Alignment/BAM

jobid_compress=$(cat $output_dir/logs/jobids.txt | grep -w "COMPRESS_FILE" | cut -d ':' -f2 )

if [ -d $output_dir/QC/trimmomatic ]; then
	input_dir=$output_dir/QC/trimmomatic
	if [ "$use_sge" = "1" ]; then
		mapping_args="${SGE_ARGS} -pe openmp $threads -hold_jid ${jobid_compress}"
	fi
else
  	input_dir=$output_dir/raw
fi

if [ $trimming == "YES" ]; then

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

if [ $mapping == "YES" ]; then
 
	if [ "$use_sge" = "1" ]; then
	
	mapping_qsub=$( qsub $mapping_args -t 1-$sample_count -N $JOBNAME.MAPPING $mapping_cmd)
	jobid_mapping=$( echo $mapping_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
	
	echo -e "MAPPING:$jobid_mapping\n" >> $output_dir/logs/jobids.txt
		samtobam_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_mapping"
		samtobam=$( qsub $samtobam_args -t 1-$sample_count -N $JOBNAME.SAMTOBAM $samtobam_cmd)
		jobid_samtobam=$( echo $samtobam | cut -d ' ' -f3 | cut -d '.' -f1 )
	
	echo -e "SAMTOBAM:$jobid_samtobam\n" >> $output_dir/logs/jobids.txt
	
	else                              
		for count in `seq 1 $sample_count`
	
		do                                                                                                           	
		echo "Running mapping on sample $count"
			execute_mapping=$($mapping_cmd $count)
      			samtobam=$($samtobam_cmd $count)
     		done                                                                                                                                                                                                                              fi                                                                                                                                                                                                                                fi

picard_cmd="$SCRIPTS_DIR/picard_duplicates.sh \
	$output_dir/Alignment/BAM \
	$samples \
	$mappingArray_rg_list \
	$duplicateBamArray_list \
	$picard_path"

if [ $duplicate_filter == "YES" ]; then
	if [ "$use_sge" = "1" ]; then
		picard_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_samtobam"
   		picard=$( qsub $picard_args -t 1-$sample_count -N $JOBNAME.PICARD $picard_cmd)
		jobid_picard=$( echo $picard | cut -d ' ' -f3 | cut -d '.' -f1 )
		
		echo -e "PICARD:$jobid_picard\n" >> $output_dir/logs/jobids.txt
  	else
		for count in `seq 1 $sample_count`
		do
		echo "Running mapping on sample $count"
			picard=$( $picard_cmd $count)

		done
	fi
fi

## FASTQC                                                                                                                                                                      
bamutil_preduplicates="$SCRIPTS_DIR/bamutil.sh \
	$output_dir/Alignment/BAM \
	$samples \
	$mappingArray_rg_list \
	$bamstatArray_pre_list \
	$exome_enrichement"

bamutil_postduplicates="$SCRIPTS_DIR/bamutil.sh \
	$output_dir/Alignment/BAM \
	$samples \
	$duplicateBamArray_list \
	$bamstatArray_post_list \
	$exome_enrichement"  
        
if [ "$use_sge" = "1" ]; then                                                                                                                                                  
	bamutil_pre_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_samtobam"
	bamutil_pre=$( qsub $bamutil_pre_args -t 1-$sample_count -N $JOBNAME.BAMUTIL_PRE $bamutil_preduplicates)                                                                                
        jobid_bamutil_pre=$( echo $bamutil_pre | cut -d ' ' -f3 | cut -d '.' -f1 )                                                                                                   
        echo -e "BAMUTIL_PRE:$jobid_bamutil_pre\n" >> $output_dir/logs/jobids.txt
 	
	if [ $duplicate_filter == "YES" ]; then                                                                                                                                            
 		bamutil_post_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_picard"                                                                                                                 
 		bamutil_post=$( qsub $bamutil_post_args -t 1-$sample_count -N $JOBNAME.BAMUTIL_POST $bamutil_postduplicates)                                                                       
 		jobid_bamutil_post=$( echo $bamutil_post | cut -d ' ' -f3 | cut -d '.' -f1 )
 		echo -e "BAMUTIL_POST:$jobid_bamutil_post\n" >> $output_dir/logs/jobids.txt
 	fi                                                                                                                                                                         
else                                                                                                                                                                           
	for count in `seq 1 $sample_count`
	do                                                                                                                   echo "Running fastqc on sample $count"
	bamutil_pre=$( $bamutil_preduplicates $count)                                                      
 		
		if [ $duplicate_filter == "YES" ]; then
			bamutil_post=$( $bamutil_postduplicates $count)
		fi
	done
fi
