uthor:A.Hernandez
#help
#Usage: cfsan_snp_matrix.sh 

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
consensus_fasta_list=$3
consensus_preserved_fasta_list=$4
snpma_fasta=$5
snpma_preser_fasta=$6


if [ "$use_sge" = "1" ]; then                                                                                                                                                                                               
   	sample_count=$sge_task_id                                                                      
else                                                                                                        
   	sample_count=$7                                                                                   
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
consensus_fasta_file=$( echo $consensus_fasta_list | tr ":" "\n" | head -$sample_count | tail -1)
consensus_preserved_fasta_file=$( echo $consensus_preserved_fasta_list | tr ":" "\n" | head -$sample_count | tail -1)

cfsan_snp_pipeline snp_matrix -c $dir/$sample/$consensus_fasta_file -o $dir/$snpma_fasta $dir/sampleDirectories.txt.OrigVCF.filtered
cfsan_snp_pipeline snp_matrix -c $dir/$samples/$consensus_preserved_fasta_file -o $dir/$snpma_preser_fasta $dir/sampleDirectories.txt.PresVCF.filtered
