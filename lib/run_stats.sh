#!/bin/bash
## Author A.Hernandez
## Usage: run_stats.sh ....

CONFIG_FILE=$1

#Execute processing_config.sh
if [ -z $SCRIPTS_DIR ]; then
        SCRIPTS_DIR=$( cat $CONFIG_FILE | grep -w 'SCRIPTS_DIR' | cut -d '=' -f2 )
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

else
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"
fi


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

##create directories
mkdir -p $output_dir/stats

jobid_samtobam=$(cat $output_dir/logs/jobids.txt | grep -w "SAMTOBAM" | cut -d ':' -f2 )
jobid_picard=$(cat $output_dir/logs/jobids.txt | grep -w "PICARD" | cut -d ':' -f2 )
jobid_msa_allsnp=$(cat $output_dir/logs/jobids.txt | grep -w "TSV_TO_MSAallsnp" | cut -d ':' -f2 )
jobid_msa_fil=$(cat $output_dir/logs/jobids.txt | grep -w "TSV_TO_MSAfil" | cut -d ':' -f2 )


if [ $stats == "YES" ]; then
	
	coverage_cmd="$SCRIPTS_DIR/bedtool_coverage.sh \
		$output_dir/Alignment/BAM \
		$output_dir/stats \
		$samples \
		$duplicateBamArray_list \
		$coverage_list \
		$coverage_graph_list \
		$ref_path"

	r_script_cmd="$SCRIPTS_DIR/graphs_coverage.R \
		$output_dir/stats"

	distances_allsnp_cmd="$SCRIPTS_DIR/distances.R \
		$output_dir/variant_calling/variants_gatk/variants \
                $output_dir/stats \
                $msa_allsnp_file \
		$dist_allsnp \
                $dist_pair_allsnp"
	
	distances_fil_cmd="$SCRIPTS_DIR/distances.R \
                $output_dir/variant_calling/variants_gatk/variants \
                $output_dir/stats \
                $msa_fil_file \
		$dist_fil \
		$dist_pair_fil"

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


fi

if [ "$use_sge" = "1" ]; then
	bedtool_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid_picard"
	bedtool_qsub=$( qsub $bedtool_args -t 1-$sample_count -N $JOBNAME.BEDTOOL $coverage_cmd)
	jobid_coverage=$( echo $bedtool_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "COVERAGE:$jobid_coverage\n" >> $output_dir/logs/jobids.txt	

	graph_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid_coverage"
	graph_qsub=$( qsub $graph_args -N $JOBNAME.GRAPH $r_script_cmd)
	jobid_graph=$( echo $graph_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "GRAPH:$jobid_graph\n" >> $output_dir/logs/jobids.txt
	
	distances_allsnp_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid_msa_allsnp"
        distances_allsnp_qsub=$( qsub $distances_allsnp_args -N $JOBNAME.DISTANCEallsnp $distances_allsnp_cmd)
        jobid_distances_allsnp=$( echo $distances_allsnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "DISTANCESallsnp:$jobid_distances_allsnp\n" >> $output_dir/logs/jobids.txt

	distances_fil_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid_msa_fil"
        distances_fil_qsub=$( qsub $distances_fil_args -N $JOBNAME.DISTANCEfil $distances_fil_cmd)
        jobid_distances_fil=$( echo $distances_fil_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "DISTANCESfil:$jobid_distances_fil\n" >> $output_dir/logs/jobids.txt
	
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
	do
		echo "Running stats on sample $count"
		stats=$( $coverage_cmd $count)

		echo "Running fastqc on sample $count"
        	bamutil_pre=$( $bamutil_preduplicates $count)
	
		if [ $duplicate_filter == "YES" ]; then
                        bamutil_post=$( $bamutil_postduplicates $count)
                fi

	done
	$r_script_cmd
	$distances_allsnp_cmd
	$distances_fil_cmd

fi
