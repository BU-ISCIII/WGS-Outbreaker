#!/bin/bash
## Author: A. Hernandez
# Help
# usage: raxml_inference.sh ....
#


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#VARIABLES


dir=$1
output_dir=$2
snp_msa=$3
model=$4

#if [ "$use_sge" = "1" ]; then

#	run_inference=$( qsub ${SGE_ARGS} -N $JOBNAME.RAXMLinfe -hold_jid $jobid -pe orte 100 mpirun raxmlHPC-MPI-AVX -m $model -V  -w $output_dir -n RAXML_TREE_INFERENCE -p 12345 -s $dir/$snp_msa -N 100)
	

#	jobid_raxmlinfe=$(echo $run_inference | cut -d ' ' -f3 | cut -d '.' -f1 )
	

#	run_boot=$( qsub ${SGE_ARGS} -N $JOBNAME.RAXMLboot -hold_jid $jobid_raxmlinfe -pe orte 100 mpirun raxmlHPC-MPI-AVXraxmlHPC-MPI-AVX -m $model -V -b 12345 -w $output_dir -n RAXML_TREE_BOOTSTRAP -p 12345 -s $dir/$snp_msa -N 100)

#	jobid_raxmlboot=$(echo $run_boot | cut -d ' ' -f3 | cut -d '.' -f1 )


#	 run_annot=$( qsub ${SGE_ARGS} -N $JOBNAME.RAXMLannot -hold_jid $jobid_raxmlboot -pe orte 100 mpirun raxmlHPC-AVX -f b -p 12345  -m $model  -w $output_dir -t $output_dir/RAxML_bestTree.RAXML_TREE_INFERENCE -z RAxML_bootstrap.RAXML_TREE_BOOTSTRAP -n RAXML_TREE_ANNOT)

#	 jobid_raxmlannot=$(echo $run_annot | cut -d ' ' -f3 | cut -d '.' -f1 )

	
	#jobid_raxmlinfe_cfsan_allsnp=$(echo $raxmlinfe_cfsan_allsnp_exe | cut -d ' ' -f3 | cut -d '.' -f1 )

	#echo -e "RAXMLinfe_CFSANallsnp:$jobid_raxmlinfe_cfsan_allsnp\n" >> $output_dir/logs/jobids.txt


#else

raxmlHPC-MPI-AVX -m $model -V -w $output_dir -n RAXML_TREE_INFERENCE -p 12345 -s $dir/$snp_msa -N 100

#fi


