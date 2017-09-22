uthor:A.Hernandez
#help
#Usage: cfsan_call_consensus.sh 

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

#VARIABLES

dir=$1
samples=$2
snp_list=$3
snp_list_preser=$4
consensus_fasta_list=$5
consensus_vcf_list=$6
consensus_preserved_fasta_list=$7
consensus_preserved_vcf_list=$8

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                               
   	sample_count=$sge_task_id                                                                      
else                                                                                                        
   	sample_count=$9                                                                                   
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
consensus_fasta_file=$( echo $consensus_fasta_list | tr ":" "\n" | head -$sample_count | tail -1)
consensus_vcf_file=$( echo $consensus_vcf_list | tr ":" "\n" | head -$sample_count | tail -1)
consensus_preserved_fasta_file=$( echo $consensus_preserved_fasta_list= | tr ":" "\n" | head -$sample_count | tail -1)
consensus_preserved_vcf_file=$( echo $consensus_preserved_vcf_list= | tr ":" "\n" | head -$sample_count | tail -1)

cfsan_snp_pipeline call_consensus -l $dir/$snp_list --vcfFileName $consensus_vcf_file -o $dir/$sample/$consensus_fasta_file $dir/$sample/reads.all.pileup

cfsan_snp_pipeline call_consensus -l $dir/$snp_list --vcfFileName $consensus_preserved_vcf_file -o $dir/$sample/$consensus_preserved_fasta_file $dir/$sample/reads.all.pileup

