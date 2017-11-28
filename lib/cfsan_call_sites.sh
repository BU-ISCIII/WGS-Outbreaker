#!/bin/bash
## Author:A.Hernandez
## version v2.0
## Usage: cfsan_call_sites.sh 

if [ $# -eq 0 ]; then
	echo -e "\nScript to tun cfsan call_sites\n"
	echo -e "Usage: cfsan_call_sites.sh input_dir samples_list unsortedBam_list VarScan_qual VarScan_frec samtoolsQ reference_path"
	exit
fi

# Test whether the script is being executed with sge or not.
if [ -z $SGE_TASK_ID ]; then
        use_sge=0
else
        use_sge=1
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#VARIABLES

dir=$1
samples=$2
unsorted_bam_list=$3
VarScan_qual=$4
VarScan_frec=$5
samtoolsQ=$6
cfsan_ref_path=$7

if [ "$use_sge" = "1" ]; then
	sample_count=$SGE_TASK_ID
else
	sample_count=$8                                                                             
fi

export SamtoolsMpileup_ExtraParams="-q 0 -Q $samtoolsQ"
export VarscanMpileup2snp_ExtraParams="--min-avg-qual $VarScan_qual --min-var-freq $VarScan_frec"
export VarscanJvm_ExtraParams="$JAVA_RAM"
export PicardJvm_ExtraParams="$JAVA_RAM"

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
unsorted_bam_file=$( echo $unsorted_bam_list | tr ":" "\n" | head -$sample_count | tail -1)

cfsan_snp_pipeline call_sites $cfsan_ref_path $dir/$sample
