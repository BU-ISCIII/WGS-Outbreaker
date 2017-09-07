#!/bin/bash
## Author: A.Hernandez
## Usage: processing_config.sh config_file

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


CONFIG_FILE=/home/ahernandez/outbreakWGS/config.file

# RunInfo
email=$( cat $CONFIG_FILE | grep -w 'MAIL' | cut -d '=' -f2 )

date_run=$( cat $CONFIG_FILE | grep -w 'DATE_RUN' | cut -d '=' -f2 )
platform=$( cat $CONFIG_FILE | grep -w 'PLATFORM' | cut -d '=' -f2 )
model=$( cat $CONFIG_FILE | grep -w 'MODEL' | cut -d '=' -f2 )
library=$( cat $CONFIG_FILE | grep -w 'LIBRARY' | cut -d '=' -f2 )
sequencing_center=$( cat $CONFIG_FILE | grep -w 'SEQUENCING_CENTER' | cut -d '=' -f2 )
run_platform=$( cat $CONFIG_FILE | grep -w 'RUN_PLATFORM' | cut -d '=' -f2 )

use_sge=$( cat $CONFIG_FILE | grep -w 'USE_SGE' | cut -d '=' -f2 )
samples=$( cat $CONFIG_FILE | grep -w 'SAMPLES' | cut -d '=' -f2 )

input_dir=$( cat $CONFIG_FILE | grep -w 'INPUT_DIR' | cut -d '=' -f2 )
output_dir=$( cat $CONFIG_FILE | grep -w 'OUTPUT_DIR' | cut -d '=' -f2 )
threads=$( cat $CONFIG_FILE | grep -w 'THREADS' | cut -d '=' -f2 )

# Pipeline steps
trimming=$( cat $CONFIG_FILE | grep -w 'TRIMMING' | cut -d '=' -f2 )
mapping=$( cat $CONFIG_FILE | grep -w 'MAPPING' | cut -d '=' -f2 )
duplicate_filter=$( cat $CONFIG_FILE | grep -w 'DUPLICATE_FILTER' | cut -d '=' -f2 )
variant_calling=$( cat $CONFIG_FILE | grep -w 'VARIANT_CALLING' | cut -d '=' -f2 )
kmerfinder=$( cat $CONFIG_FILE | grep -w 'KMERFINDER' | cut -d '=' -f2 )
srst2=$( cat $CONFIG_FILE | grep -w 'SRST2' | cut -d '=' -f2 )

# REFERENCES
exome_enrichement=$( cat $CONFIG_FILE | grep -w 'EXOME_ENRICHMENT' | cut -d '=' -f2 )
ref_path=$( cat $CONFIG_FILE | grep -w 'GENOME_REF' | cut -d '=' -f2 )
know_snps=$( cat $CONFIG_FILE | grep -w 'KNOWN_SNPS' | cut -d '=' -f2 )
know_indels=$( cat $CONFIG_FILE | grep -w 'KNOWN_INDELS' | cut -d '=' -f2 )
bact_db_path=$( cat $CONFIG_FILE | grep -w 'BACT_DB_PATH' | cut -d '=' -f2 )
srst2_db_path_argannot=$( cat $CONFIG_FILE | grep -w 'SRST2_DB_PATH_ARGannot' | cut -d '=' -f2)
srst2_db_path_plasmidfinder=$( cat $CONFIG_FILE | grep -w 'SRST2_DB_PATH_PlasmidFinder' | cut -d '=' -f2)
srst2_db_path_mlst_db=$( cat $CONFIG_FILE | grep -w 'SRST2_DB_PATH_mlst_db' | cut -d '=' -f2)
srst2_db_path_mlst_definitions=$( cat $CONFIG_FILE | grep -w 'SRST2_DB_PATH_mlst_definitions' | cut -d '=' -f2)


# Arguments
trimmomatic_version=$( cat $CONFIG_FILE | grep -w 'trimmomatic_version' | cut -d '=' -f2 )
trimmomatic_path=$( cat $CONFIG_FILE | grep -w 'TRIMMOMATIC_PATH' | cut -d '=' -f2 )
trim_args=$( cat $CONFIG_FILE | grep -w 'TRIM_ARGS' | cut -d '=' -f2 )
picard_path=$( cat $CONFIG_FILE | grep -w 'PICARD_PATH' | cut -d '=' -f2 )
gatk_path=$( cat $CONFIG_FILE | grep -w 'GATK_PATH' | cut -d '=' -f2 )
kmerfinder_path=$( cat $CONFIG_FILE | grep -w 'KMERFINDER_PATH' | cut -d '=' -f2 )
#SRST2_DELIMITER=$( cat $CONFIG_FILE | grep -w 'SRST2_DELIMITER' | cut -d '=' -f2 )


## Extract fastq and names for samples.
sample_count=$(echo $samples | tr ":" "\n" | wc -l)

# Collect fastq files
mkdir -p $output_dir/raw
i=1
for sample in $(echo $samples | tr ":" " ")
do
        mkdir -p $output_dir/raw/$sample
        ## Raw reads Array
        fastqArray_R1[$i]=$( grep -w "^${sample}" $CONFIG_FILE | cut -d '=' -f2 | cut -d$'\t' -f1 )
        fastqArray_R2[$i]=$( grep -w "^${sample}" $CONFIG_FILE | cut -d '=' -f2 | cut -d$'\t' -f2 )

        # Create raw folder with raw fastq links
        ln -fs $input_dir/${fastqArray_R1[$i]} $output_dir/raw/$sample/${fastqArray_R1[$i]}
        ln -fs $input_dir/${fastqArray_R2[$i]} $output_dir/raw/$sample/${fastqArray_R2[$i]}

        # Create trimming names
        trimmedFastqArray_paired_R1[$i]=$sample.trimmed_R1.fastq
        trimmedFastqArray_paired_R2[$i]=$sample.trimmed_R2.fastq
        trimmedFastqArray_unpaired_R1[$i]=$sample.trimmed_unpaired_R1.fastq
        trimmedFastqArray_unpaired_R2[$i]=$sample.trimmed_unpaired_R2.fastq

        #create compress name
        compress_paired_R1[$i]=$sample.trimmed_R1.fastq.gz
        compress_paired_R2[$i]=$sample.trimmed_R2.fastq.gz
        compress_unpaired_R1[$i]=$sample.trimmed_unpaired__R1.fastq.gz
        compress_unpaired_R2[$i]=$sample.trimmed_unpaired_R2.fastq.gz

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

        #Create concat output names
        concatFastq[$i]=$sample.concat.fastq

        #create kmerfinder output names
        kmerfinderST[$i]=$sample.kmerfinder.txt

        #create srst2 output names
        resistance[$i]=$sample.resistance
        plasmid[$i]=$sample.plasmid
	mlst[$i]=$sample.mlst

        let i=i+1
done

vcfArray_list=all_samples.vcf
vcfsnpsArray_list=all_samples_snps.vcf
vcfsnpsfilArray_list=all_samples_snps_fil.vcf
vcfindelsArray_list=all_samples_indels.vcf
vcfindelsfilArray_list=all_samples_indels_fil.vcf
vcffilArray_list=all_samples_fil.vcf

fastq_R1_list=$( echo ${fastqArray_R1[@]} | tr " " ":" )
fastq_R2_list=$( echo ${fastqArray_R2[@]} | tr " " ":" )
trimmedFastqArray_paired_R1_list=$( echo ${trimmedFastqArray_paired_R1[@]} | tr " " ":" )
trimmedFastqArray_paired_R2_list=$( echo ${trimmedFastqArray_paired_R2[@]} | tr " " ":" )
trimmedFastqArray_unpaired_R1_list=$( echo ${trimmedFastqArray_unpaired_R1[@]} | tr " " ":" )
trimmedFastqArray_unpaired_R2_list=$( echo ${trimmedFastqArray_unpaired_R2[@]} | tr " " ":" )
compress_paired_R1_list=$( echo ${compress_paired_R1[@]} | tr " " ":" )
compress_paired_R2_list=$( echo ${compress_paired_R2[@]} | tr " " ":" )
compress_unpaired_R1_list=$( echo ${compress_unpaired_R1[@]} | tr " " ":" )
compress_unpaired_R2_list=$( echo ${compress_unpaired_R2[@]} | tr " " ":" )
mappingArray_sam_list=$( echo ${mappingArray_sam[@]} | tr " " ":" )
mappingArray_bam_list=$( echo ${mappingArray_bam[@]} | tr " " ":" )
mappingArray_sorted_list=$( echo ${mappingArray_sorted[@]} | tr " " ":" )
mappingArray_rg_list=$( echo ${mappingArray_rg[@]} | tr " " ":" )
bamstatArray_pre_list=$( echo ${bamstatArray_pre[@]} | tr " " ":" )
bamstatArray_post_list=$( echo ${bamstatArray_post[@]} | tr " " ":" )
duplicateBamArray_list=$( echo ${duplicateBamArray[@]} | tr " " ":" )
recalibratedBamArray_list=$( echo ${recalibratedBamArray[@]} | tr " " ":" )
realignedBamArray_list=$( echo ${realignedBamArray[@]} | tr " " ":" )
concatFastq_list=$( echo ${concatFastq[@]} | tr " " ":")
kmerfinderST_list=$( echo ${kmerfinderST[@]} | tr " " ":")
resistance_list=$( echo ${resistance[@]} | tr " " ":")
plasmid_list=$( echo ${plasmid[@]} | tr " " ":")
mlst_list=$( echo ${mlst[@]} | tr " " ":")
