#!/bin/bash
##Auhtor: A.Hernandez

#help
#Usage: cfsan_snp_list.sh 

# Test whether the script is being executed with sge or not.
if [ -z $sge_task_id ]; then
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

#VARAIBLES

var_vcf=$1
var_preserv_vcf=$2
dir=$3


cfsan_snp_pipeline merge_sites -n $var_vcf -o $dir/snplist.txt $dir/sampleDirectories.txt $dir/sampleDirectories.txt.OrigVCF.filtered
cfsan_snp_pipeline merge_sites -n $var_preserv_vcf -o $dir/snplist_preserved.txt $dir/sampleDirectories.txt $dir/sampleDirectories.txt.PresVCF.filtered
