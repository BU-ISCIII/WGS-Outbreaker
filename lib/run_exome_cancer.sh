#!/bin/bash        
## Author S. Monzon
## version v2.0    
## usage: run_exome.sh <config.file>

###############
## VARIABLES ##
###############
set -e
set -u
set -x

# Configuration
CONFIG_FILE=$1

#Global VARIABLES
export JOBNAME="exome_pipeline_v2.0"
export SCRIPTS_DIR=$( cat $CONFIG_FILE | grep -w 'SCRIPTS_DIR' | cut -d '=' -f2 )
export TEMP=$( cat $CONFIG_FILE | grep -w 'TEMP_DIR' | cut -d '=' -f2 )
export JAVA_RAM=$( cat $CONFIG_FILE | grep -w 'JAVA_RAM' | cut -d '=' -f2 )  

# RunInfo
EMAIL=$( cat $CONFIG_FILE | grep -w 'MAIL' | cut -d '=' -f2 )

DATE_RUN=$( cat $CONFIG_FILE | grep -w 'DATE_RUN' | cut -d '=' -f2 ) 
PLATFORM=$( cat $CONFIG_FILE | grep -w 'PLATFORM' | cut -d '=' -f2 ) 
MODEL=$( cat $CONFIG_FILE | grep -w 'MODEL' | cut -d '=' -f2 ) 
LIBRARY=$( cat $CONFIG_FILE | grep -w 'LIBRARY' | cut -d '=' -f2 ) 
SEQUENCING_CENTER=$( cat $CONFIG_FILE | grep -w 'SEQUENCING_CENTER' | cut -d '=' -f2 ) 
RUN_PLATFORM=$( cat $CONFIG_FILE | grep -w 'RUN_PLATFORM' | cut -d '=' -f2 )  

USE_SGE=$( cat $CONFIG_FILE | grep -w 'USE_SGE' | cut -d '=' -f2 ) 
SAMPLES=$( cat $CONFIG_FILE | grep -w 'SAMPLES' | cut -d '=' -f2 )
CONTROL=$( cat $CONFIG_FILE | grep -w 'CONTROL' | cut -d '=' -f2 ) 
CASE=$( cat $CONFIG_FILE | grep -w 'CASE' | cut -d '=' -f2 )   

INPUT_DIR=$( cat $CONFIG_FILE | grep -w 'INPUT_DIR' | cut -d '=' -f2 ) 
OUTPUT_DIR=$( cat $CONFIG_FILE | grep -w 'OUTPUT_DIR' | cut -d '=' -f2 ) 
THREADS=$( cat $CONFIG_FILE | grep -w 'THREADS' | cut -d '=' -f2 ) 

# Pipeline steps
TRIMMING=$( cat $CONFIG_FILE | grep -w 'TRIMMING' | cut -d '=' -f2 )  
MAPPING=$( cat $CONFIG_FILE | grep -w 'MAPPING' | cut -d '=' -f2 )   
DUPLICATE_FILTER=$( cat $CONFIG_FILE | grep -w 'DUPLICATE_FILTER' | cut -d '=' -f2 )    
VARIANT_CALLING=$( cat $CONFIG_FILE | grep -w 'VARIANT_CALLING' | cut -d '=' -f2 )     
VARIANT_CALLER=$( cat $CONFIG_FILE | grep -w 'VARIANT_CALLER' | cut -d '=' -f2 )      

# REFERENCES
REF_PATH=$( cat $CONFIG_FILE | grep -w 'GENOME_REF' | cut -d '=' -f2 ) 
EXOME_ENRICHMENT=$( cat $CONFIG_FILE | grep -w 'EXOME_ENRICHMENT' | cut -d '=' -f2 )

# Arguments
trimmomatic_version=$( cat $CONFIG_FILE | grep -w 'trimmomatic_version' | cut -d '=' -f2 ) 
TRIMMOMATIC_PATH=$( cat $CONFIG_FILE | grep -w 'TRIMMOMATIC_PATH' | cut -d '=' -f2 ) 
TRIM_ARGS=$( cat $CONFIG_FILE | grep -w 'TRIM_ARGS' | cut -d '=' -f2 )
PICARD_PATH=$( cat $CONFIG_FILE | grep -w 'PICARD_PATH' | cut -d '=' -f2 ) 
STRELKA_CONFIG=$( cat $CONFIG_FILE | grep -w 'STRELKA_CONFIG' | cut -d '=' -f2 )
GATK_PATH=$( cat $CONFIG_FILE | grep -w 'GATK_PATH' | cut -d '=' -f2 )  

## Extract fastq and names for samples.
sample_count=$(echo $SAMPLES | tr ":" "\n" | wc -l)
pair_count=$( echo $CASE | tr ":" "\n" | wc -l)

## SGE args
if [ "$USE_SGE"="1" ]; then
  mkdir -p $OUTPUT_DIR/logs
  export SGE_ARGS="-V -j y -b y -wd $OUTPUT_DIR/logs -m a -M $EMAIL"
fi

# Collect fastq files
mkdir -p $OUTPUT_DIR/raw
i=1
for sample in $(echo $SAMPLES | tr ":" " ")
do
	mkdir -p $OUTPUT_DIR/raw/$sample
	## Raw reads Array
	fastqArray_R1[$i]=$( grep -w "^${sample}" $CONFIG_FILE | cut -d '=' -f2 | cut -d$'\t' -f1 )
	fastqArray_R2[$i]=$( grep -w "^${sample}" $CONFIG_FILE | cut -d '=' -f2 | cut -d$'\t' -f2 )
	
	# Create raw folder with raw fastq links
	ln -fs $INPUT_DIR/${fastqArray_R1[$i]} $OUTPUT_DIR/raw/$sample/${fastqArray_R1[$i]} 
 	ln -fs $INPUT_DIR/${fastqArray_R2[$i]} $OUTPUT_DIR/raw/$sample/${fastqArray_R2[$i]}  
	
    # Create trimming names
    trimmedFastqArray_paired_R1[$i]=$sample.trimmed_R1.fastq            
    trimmedFastqArray_paired_R2[$i]=$sample.trimmed_R2.fastq            
    trimmedFastqArray_unpaired_R1[$i]=$sample.trimmed_unpaired_R1.fastq 
    trimmedFastqArray_unpaired_R2[$i]=$sample.trimmed_unpaired_R2.fastq 

	# Create mapping names
	mappingArray_sam[$i]=$sample.align.sam
	mappingArray_bam[$i]=$sample.align.bam
	mappingArray_sorted[$i]=$sample.align.sorted.bam
    mappingArray_rg[$i]=$sample.align.sorted.rg.bam

	# Create bamstat names
	bamstatArray_pre[$i]=$sample.pre.bamstat.txt
	bamstatArray_post[$i]=$sample.post.bamstat.txt

	# Create duplicate bam names
	duplicateBamArray[$i]=$sample.woduplicates.bam

	# Create gatk output names
	recalibratedBamArray[$i]=$sample.recalibrated.bam
	realignedBamArray[$i]=$sample.realigned.bam

	let i=i+1
done

if [ $VARIANT_CALLER == "STRELKA" ];then
	for tumor_sample in $(echo $CASE)
	do
		# Vcf names
   	 	vcfArray[$i]=$tumor_sample.vcf
		let i=i+1
	done
 	vcfArray_list=$( echo ${vcfArray[@]} | tr " " ":" ) 
fi

if [ $VARIANT_CALLER == "GATK_HAPLOID" ];then
   	vcfArray_list=all_samples.vcf
   	vcfsnpsArray_list=all_samples_snps.vcf
   	vcfsnpsfilArray_list=all_samples_snps_fil.vcf
   	vcfindelsArray_list=all_samples_indels.vcf
   	vcfindelsfilArray_list=all_samples_indels_fil.vcf
fi	

fastq_R1_list=$( echo ${fastqArray_R1[@]} | tr " " ":" )
fastq_R2_list=$( echo ${fastqArray_R2[@]} | tr " " ":" )
trimmedFastqArray_paired_R1_list=$( echo ${trimmedFastqArray_paired_R1[@]} | tr " " ":" )       
trimmedFastqArray_paired_R2_list=$( echo ${trimmedFastqArray_paired_R2[@]} | tr " " ":" )       
trimmedFastqArray_unpaired_R1_list=$( echo ${trimmedFastqArray_unpaired_R1[@]} | tr " " ":" )   
trimmedFastqArray_unpaired_R2_list=$( echo ${trimmedFastqArray_unpaired_R2[@]} | tr " " ":" )   
mappingArray_sam_list=$( echo ${mappingArray_sam[@]} | tr " " ":" )
mappingArray_bam_list=$( echo ${mappingArray_bam[@]} | tr " " ":" )
mappingArray_sorted_list=$( echo ${mappingArray_sorted[@]} | tr " " ":" )
mappingArray_rg_list=$( echo ${mappingArray_rg[@]} | tr " " ":" ) 
bamstatArray_pre_list=$( echo ${bamstatArray_pre[@]} | tr " " ":" )
bamstatArray_post_list=$( echo ${bamstatArray_post[@]} | tr " " ":" ) 
duplicateBamArray_list=$( echo ${duplicateBamArray[@]} | tr " " ":" ) 
recalibratedBamArray_list=$( echo ${recalibratedBamArray[@]} | tr " " ":" ) 
realignedBamArray_list=$( echo ${realignedBamArray[@]} | tr " " ":" )   


# Execute preprocessing
$SCRIPTS_DIR/run_preprocessing.sh $TRIMMING $USE_SGE $INPUT_DIR $OUTPUT_DIR $THREADS $fastq_R1_list $fastq_R2_list $sample_count $SAMPLES $TRIM_ARGS $trimmomatic_version $TRIMMOMATIC_PATH $trimmedFastqArray_paired_R1_list $trimmedFastqArray_paired_R2_list $trimmedFastqArray_unpaired_R1_list $trimmedFastqArray_unpaired_R2_list

# Execute MAPPING
if [ $TRIMMING == "YES" ]; then
    $SCRIPTS_DIR/run_mapping.sh $MAPPING $USE_SGE $OUTPUT_DIR $REF_PATH $THREADS $trimmedFastqArray_paired_R1_list $trimmedFastqArray_paired_R2_list $sample_count $SAMPLES $mappingArray_sam_list $mappingArray_bam_list $mappingArray_sorted_list $bamstatArray_pre_list $bamstatArray_post_list $duplicateBamArray_list $DUPLICATE_FILTER $PICARD_PATH $EXOME_ENRICHMENT $PLATFORM $MODEL $DATE_RUN $LIBRARY $SEQUENCING_CENTER $RUN_PLATFORM $mappingArray_rg_list
else
	$SCRIPTS_DIR/run_mapping.sh $MAPPING $USE_SGE $OUTPUT_DIR $REF_PATH $THREADS $fastq_R1_list $fastq_R2_list $sample_count $SAMPLES $mappingArray_sam_list $mappingArray_bam_list $mappingArray_sorted_list $bamstatArray_pre_list $bamstatArray_post_list $duplicateBamArray_list $DUPLICATE_FILTER $PICARD_PATH $EXOME_ENRICHMENT $PLATFORM $MODEL $DATE_RUN $LIBRARY $SEQUENCING_CENTER $RUN_PLATFORM $mappingArray_rg_list 
fi

# Execute variant Calling
if [ $DUPLICATE_FILTER == "YES" ]; then
	$SCRIPTS_DIR/run_variantCalling_cancer.sh $USE_SGE $VARIANT_CALLING $VARIANT_CALLER $DUPLICATE_FILTER $OUTPUT_DIR $REF_PATH $THREADS $SAMPLES $duplicateBamArray_list $EXOME_ENRICHMENT $CONTROL $CASE $vcfArray_list $STRELKA_CONFIG $pair_count
else
    $SCRIPTS_DIR/run_variantCalling_cancer.sh $USE_SGE $VARIANT_CALLING $VARIANT_CALLER $DUPLICATE_FILTER $OUTPUT_DIR $REF_PATH $THREADS $SAMPLES $mappingArray_sorted_list $EXOME_ENRICHMENT $CONTROL $CASE $vcfArray_list $STRELKA_CONFIG $pair_count
fi
