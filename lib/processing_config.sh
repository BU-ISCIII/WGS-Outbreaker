#!/bin/bash
## Author: A.Hernandez
## version v2.0

if [ $# -eq 0  ]; then
	echo -e "\nRead config file and set all variables\n"
        echo -e "Usage: processing_config.sh <config.file>"
        exit
fi


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#Import config file
SOURCE_ARG=$1
CONFIG_FILE=$( echo $SOURCE_ARG | cut -d '-' -f3)

#Global VARIABLES

export TEMP=$( cat $CONFIG_FILE | grep -w 'TEMP_DIR' | cut -d '=' -f2 )
export JAVA_RAM=$( cat $CONFIG_FILE | grep -w 'JAVA_RAM' | cut -d '=' -f2 )

# RunInfo
email=$( cat $CONFIG_FILE | grep -w 'MAIL' | cut -d '=' -f2 )
JOBNAME="WGS-Outbreaker_v2.0"
date_run=$( cat $CONFIG_FILE | grep -w 'DATE_RUN' | cut -d '=' -f2 )
platform=$( cat $CONFIG_FILE | grep -w 'PLATFORM' | cut -d '=' -f2 )
model=$( cat $CONFIG_FILE | grep -w 'MODEL' | cut -d '=' -f2 )
library=$( cat $CONFIG_FILE | grep -w 'LIBRARY' | cut -d '=' -f2 )
sequencing_center=$( cat $CONFIG_FILE | grep -w 'SEQUENCING_CENTER' | cut -d '=' -f2 )
run_platform=$( cat $CONFIG_FILE | grep -w 'RUN_PLATFORM' | cut -d '=' -f2 )

#HPC variables
use_sge=$( cat $CONFIG_FILE | grep -w 'USE_SGE' | cut -d '=' -f2 )
vmem=$( cat $CONFIG_FILE | grep -w 'H_VMEM' | cut -d '=' -f2)
threads=$( cat $CONFIG_FILE | grep -w 'THREADS' | cut -d '=' -f2 )
queue=$( cat $CONFIG_FILE | grep -w 'QUEUE' | cut -d '=' -f2 )

#Working directories
samples=$( cat $CONFIG_FILE | grep -w 'SAMPLES' | cut -d '=' -f2 )
input_dir=$( cat $CONFIG_FILE | grep -w 'INPUT_DIR' | cut -d '=' -f2 )
output_dir=$( cat $CONFIG_FILE | grep -w 'OUTPUT_DIR' | cut -d '=' -f2 )

# Pipeline steps
trimming=$( cat $CONFIG_FILE | grep -w 'TRIMMING' | cut -d '=' -f2 )
check_references=$( cat $CONFIG_FILE | grep -w 'CHECK_REFERENCES' | cut -d '=' -f2)
mapping=$( cat $CONFIG_FILE | grep -w 'MAPPING' | cut -d '=' -f2 )
duplicate_filter=$( cat $CONFIG_FILE | grep -w 'DUPLICATE_FILTER' | cut -d '=' -f2 )
variant_calling=$( cat $CONFIG_FILE | grep -w 'VARIANT_CALLING' | cut -d '=' -f2 )
kmerfinder=$( cat $CONFIG_FILE | grep -w 'KMERFINDER' | cut -d '=' -f2 )
srst2=$( cat $CONFIG_FILE | grep -w 'SRST2' | cut -d '=' -f2 )
cfsan=$( cat $CONFIG_FILE | grep -w 'CFSAN' | cut -d '=' -f2 )
vcf_to_msa=$( cat $CONFIG_FILE | grep -w 'VCF_TO_MSA' | cut -d '=' -f2 )
raxml=$( cat $CONFIG_FILE | grep -w 'RAXML' | cut -d '=' -f2 )
stats=$( cat $CONFIG_FILE | grep -w 'STATS' | cut -d '=' -f2 )


# REFERENCES
exome_enrichement=$( cat $CONFIG_FILE | grep -w 'EXOME_ENRICHMENT' | cut -d '=' -f2 )
ref_path=$( cat $CONFIG_FILE | grep -w 'GENOME_REF' | cut -d '=' -f2 )
genome_name=$( cat $CONFIG_FILE | grep -w 'GENOME_NAME' | cut -d '=' -f2 )
know_snps=$( cat $CONFIG_FILE | grep -w 'KNOWN_SNPS' | cut -d '=' -f2 )
know_indels=$( cat $CONFIG_FILE | grep -w 'KNOWN_INDELS' | cut -d '=' -f2 )
bact_db_path=$( cat $CONFIG_FILE | grep -w 'BACT_DB_PATH' | cut -d '=' -f2 )
srst2_db_path_argannot=$( cat $CONFIG_FILE | grep -w 'SRST2_DB_PATH_ARGannot' | cut -d '=' -f2)
srst2_db_path_plasmidfinder=$( cat $CONFIG_FILE | grep -w 'SRST2_DB_PATH_PlasmidFinder' | cut -d '=' -f2)
srst2_db_path_mlst_db=$( cat $CONFIG_FILE | grep -w 'SRST2_DB_PATH_mlst_db' | cut -d '=' -f2)
srst2_db_path_mlst_definitions=$( cat $CONFIG_FILE | grep -w 'SRST2_DB_PATH_mlst_definitions' | cut -d '=' -f2)

# Arguments
#trimming
trimmomatic_version=$( cat $CONFIG_FILE | grep -w 'trimmomatic_version' | cut -d '=' -f2)
trimmomatic_path=$( cat $CONFIG_FILE | grep -w 'TRIMMOMATIC_PATH' | cut -d '=' -f2 )
trim_args=$( cat $CONFIG_FILE | grep -w 'TRIM_ARGS' | cut -d '=' -f2 )

#srst
srst2_delim=$( cat $CONFIG_FILE | grep -w 'SRST2_DELIMITER' | cut -d '=' -f2)

#cfsan
VarScan_qual=$( cat $CONFIG_FILE | grep -w 'VarScan_qual' | cut -d '=' -f2)
VarScan_frec=$( cat $CONFIG_FILE | grep -w 'VarScan_frec' | cut -d '=' -f2)
samtoolsQ=$( cat $CONFIG_FILE | grep -w 'samtoolsQ' | cut -d '=' -f2)
edge_length=$( cat $CONFIG_FILE | grep -w 'edge_length' | cut -d '=' -f2)
minBaseQual=$( cat $CONFIG_FILE | grep -w 'minBaseQual' | cut -d '=' -f2)
minConsFrec=$( cat $CONFIG_FILE | grep -w 'minConsFrec' | cut -d '=' -f2)
minConsStrdDpth=$( cat $CONFIG_FILE | grep -w 'minConsStrdDpth' | cut -d '=' -f2)
minConsStrdBias=$( cat $CONFIG_FILE | grep -w 'minConsStrdBias' | cut -d '=' -f2)

#snp filter
max_snp=$( cat $CONFIG_FILE | grep -w 'MAX_SNP' | cut -d '=' -f2)
window_size=$( cat $CONFIG_FILE | grep -w 'WINDOW_SIZE' | cut -d '=' -f2)

#Rdepth
depth=$( cat $CONFIG_FILE | grep -w 'DEPTH_COVERAGE' | cut -d '=' -f2)

# Raxml
boots=$( cat $CONFIG_FILE | grep -w 'BOOTSTRAP' | cut -d '=' -f2)
model_raxml=$( cat $CONFIG_FILE | grep -w 'MODEL_RAXML' | cut -d '=' -f2)

# paths
picard_path=$( cat $CONFIG_FILE | grep -w 'PICARD_PATH' | cut -d '=' -f2 )
gatk_path=$( cat $CONFIG_FILE | grep -w 'GATK_PATH' | cut -d '=' -f2 )
kmerfinder_path=$( cat $CONFIG_FILE | grep -w 'KMERFINDER_PATH' | cut -d '=' -f2 )


## SGE args
if [ "$use_sge" = "1" ]; then
        mkdir -p $output_dir/logs

        if [ -z "$queue" ]; then
        SGE_ARGS="-V -j y -b y -wd $output_dir/logs -m a -M $email -q all.q"
        else
        SGE_ARGS="-V -j y -b y -wd $output_dir/logs -m a -M $email $queue"
        fi
fi

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
        compress_unpaired_R1[$i]=$sample.trimmed_unpaired_R1.fastq.gz
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
        haplotypeGVCF[$i]=$sample.g.vcf
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

	#crate CFSAN output names
	align[$i]=$sample.read.sam
	sort_sam[$i]=$sample.sorted.sam
	dedup_sam[$i]=$sample.dedup.sam
	dedup_metrics[$i]=$sample.dedup.metrix.txt
	dedup_bam[$i]=$sample.dedup.bam
	dedup_bam_bai[$i]=$sample.dedup.bam.bai
	intervals[$i]=$sample.intervals
	unsorted_bam[$i]=$sample.unsorted.bam
	unsorted_bai[$i]=$sample.unsorted.bai
	var_flt[$i]=$sample.var.flv.vcf
	sort_bam[$i]=$sample.sorted.bam
	pileup[$i]=$sample.pileup
	snp_preserved[$i]=$sample.var.flt.preserved.vcf
	snp_removed[$i]=$sample.var.flt.removed.vcf
	consensus_fasta[$i]=$sample.consensus.fasta
	consensus_vcf[$i]=$sample.consensus.vcf
	consensus_preserved_fasta[$i]=$sample.consensus_preserved.fasta
        consensus_preserved_vcf[$i]=$sample.consensus_preserved.vcf
	metrics[$i]=$sample.metrics

	#create coverage stats names
	coverage[$i]=$sample.coverage.csv
	coverage_graph[$i]=$sample.coverage.graph.csv


        let i=i+1
done

#VCF FILES

#Merge GCVFs
vcfArray_list=snp_indels.vcf

#selectVariants
vcfsnpsArray_list=snp_only.vcf
vcfindelsArray_list=indels_only.vcf

#Variant filtration
vcfsnpsfilArray_list=snp_only_flags.vcf
vcfindelsfilArray_list=indels_only_flags.vcf

#VCFtool PASS and SNPCluster
vcfsnpPassCluster=snp_only_PassCluster.vcf

#VCFtool PASS snp
vcfsnpPass=snp_only_Pass.vcf

#snp and indels fil
vcffilArray_list=snp_indels_flags.vcf

#CFSAN FASTA FILES

cfsan_snpma_fasta=snpma.fasta
cfsan_snpma_fil_fasta=snpma_preserved.fasta

#VCF TO MSA FILES

tsv_filsnp_file=snp_PassCluster.tsv
msa_filsnp_file=snp_PassCluster.fasta

tsv_passnp_file=snp_Pass.tsv
msa_passnp_file=snp_Pass.fasta

#DISTANCE MATRIX FILES

dist_filsnp=distance_snp_PassCluster.txt
dist_passnp=distance_snp_Pass.txt

dist_pair_filsnp=distance_pairs_snp_PassCluster.txt
dist_pair_passnp=distance_pairs_snp_Pass.txt


#LIST FILES
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

haplotypeGVCF_list=$( echo ${haplotypeGVCF[@]} |tr " " ":")
recalibratedBamArray_list=$( echo ${recalibratedBamArray[@]} | tr " " ":" )
realignedBamArray_list=$( echo ${realignedBamArray[@]} | tr " " ":" )

concatFastq_list=$( echo ${concatFastq[@]} | tr " " ":")
kmerfinderST_list=$( echo ${kmerfinderST[@]} | tr " " ":")

resistance_list=$( echo ${resistance[@]} | tr " " ":")
plasmid_list=$( echo ${plasmid[@]} | tr " " ":")
mlst_list=$( echo ${mlst[@]} | tr " " ":")

align_list=$( echo ${align[@]} | tr " " ":")
sort_sam_list=$( echo ${sort_sam[@]} | tr " " ":")
dedup_sam_list=$( echo ${dedup_sam[@]} | tr " " ":")
dedup_metrics_list=$( echo ${dedup_metrics[@]} | tr " " ":")
dedup_bam_list=$( echo ${dedup_bam[@]} | tr " " ":")
dedup_bam_bai_list=$( echo ${dedup_bam_bai[@]} | tr " " ":")
intervals_list=$( echo ${intervals[@]} | tr " " ":")
unsorted_bam_list=$( echo ${unsorted_bam[@]} | tr " " ":")
unsorted_bai_list=$( echo ${unsorted_bai[@]} | tr " " ":")
var_flt_list=$( echo ${var_flt[@]} | tr " " ":")
sort_bam_list=$( echo ${sort_bam[@]} | tr " " ":")
pileup_list=$( echo ${pileup[@]} | tr " " ":")
snp_preserved_list=$( echo ${snp_preserved[@]} | tr " " ":")
snp_removed_list=$( echo ${snp_removed[@]} | tr " " ":")
metrics_list=$( echo ${metrics[@]} | tr " " ":")

coverage_list=$( echo ${coverage[@]} | tr " " ":")
coverage_graph_list=$( echo ${coverage_graph[@]} | tr " " ":")
