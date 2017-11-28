#!/bin/bash
## Author:A.Hernandez
## version v2.0

if [ $# -eq 0 ]; then
        echo -e "\nScrip to run cfsan call_consenssus\n"
        echo -e "Usage: cfsan_call_consensus.sh input_dir samples_list minBaseQual minConsFrec minConsStrdDpth minConsStrdBias"
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
minBaseQual=$3
minConsFrec=$4
minConsStrdDpth=$5
minConsStrdBias=$6

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                               
   	sample_count=$SGE_TASK_ID                                                                    
else                                                                                                        
   	sample_count=$7                                                                                
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)

cfsan_snp_pipeline call_consensus -l $dir/snplist.txt --minBaseQual $minBaseQual --minConsFreq $minConsFrec --minConsStrdDpth $minConsStrdDpth --minConsStrdBias $minConsStrdBias --vcfFileName consensus.vcf -o $dir/samples/$sample/consensus.fasta $dir/samples/$sample/reads.all.pileup

cfsan_snp_pipeline call_consensus -l $dir/snplist_preserved.txt --minBaseQual $minBaseQual --minConsFreq $minConsFrec --minConsStrdDpth $minConsStrdDpth --minConsStrdBias $minConsStrdBias --vcfFileName consensus_preserved.vcf -o $dir/samples/$sample/consensus_preserved.fasta -e $dir/samples/$sample/var.flt_removed.vcf $dir/samples/$sample/reads.all.pileup
