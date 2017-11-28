#!/bin/bash
## Author A.Hernandez
## version v2.0

if [ $# -eq 0  ]; then
	echo -e "\n Execute final stats\n"
        echo "Usage: run_stats.sh <config.file>"
        exit
fi

#Execute processing_config.sh

CONFIG_FILE=$1

# Check if run_outbreak_wgs.sh was execute
if [ -z $SCRIPTS_DIR ]; then
        SCRIPTS_DIR=$( cat $CONFIG_FILE | grep -w 'SCRIPTS_DIR' | cut -d '=' -f2 )
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

# Or other runner was execute
else
        source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

# Folder creation
mkdir -p $output_dir/stats

# Get jobids from samToba, duplicate filer and vcfTOmsa steps
jobid_samtobam=$(cat $output_dir/logs/jobids.txt | grep -w "SAMTOBAM" | cut -d ':' -f2 )
jobid_picard=$(cat $output_dir/logs/jobids.txt | grep -w "PICARD_DUPLICATES" | cut -d ':' -f2 )
jobid_msa_filsnp=$(cat $output_dir/logs/jobids.txt | grep -w "TSV_TO_MSAfilsnp" | cut -d ':' -f2 )
jobid_msa_passnp=$(cat $output_dir/logs/jobids.txt | grep -w "TSV_TO_MSApassnp" | cut -d ':' -f2 )

	
coverage_cmd="$SCRIPTS_DIR/bedtool_coverage.sh \
	$output_dir/Alignment/BAM \
	$output_dir/stats \
	$samples \
	$duplicateBamArray_list \
	$coverage_list \
	$coverage_graph_list \
	$ref_path"

r_script_cmd="$SCRIPTS_DIR/graphs_coverage.R \
	$output_dir/stats
	$depth"

distances_filsnp_cmd="$SCRIPTS_DIR/distances.R \
	$output_dir/variant_calling/variants_gatk/variants \
	$output_dir/stats \
	$msa_filsnp_file \
	$dist_filsnp \
	$dist_pair_filsnp"
	
distances_passnp_cmd="$SCRIPTS_DIR/distances.R \
	$output_dir/variant_calling/variants_gatk/variants \
        $output_dir/stats \
        $msa_passnp_file \
	$dist_passnp \
	$dist_pair_passnp"

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

# Execute stats in HPC
if [ "$use_sge" = "1" ]; then
	bedtool_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid_picard"
	bedtool_qsub=$( qsub $bedtool_args -t 1-$sample_count -N $JOBNAME.BEDTOOL $coverage_cmd)
	jobid_coverage=$( echo $bedtool_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "COVERAGE:$jobid_coverage\n" >> $output_dir/logs/jobids.txt	

	graph_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid_coverage"
	graph_qsub=$( qsub $graph_args -N $JOBNAME.GRAPH $r_script_cmd)
	jobid_graph=$( echo $graph_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "GRAPH:$jobid_graph\n" >> $output_dir/logs/jobids.txt
	
	distances_filsnp_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid_msa_filsnp"
        distances_filsnp_qsub=$( qsub $distances_filsnp_args -N $JOBNAME.DISTANCEfilsnp $distances_filsnp_cmd)
        jobid_distances_filsnp=$( echo $distances_filsnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "DISTANCESfilsnp:$jobid_distances_filsnp\n" >> $output_dir/logs/jobids.txt

	distances_passnp_args="${SGE_ARGS} -pe orte $threads -hold_jid $jobid_msa_passnp"
        distances_passnp_qsub=$( qsub $distances_passnp_args -N $JOBNAME.DISTANCEpassnp $distances_passnp_cmd)
        jobid_distances_passnp=$( echo $distances_passnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "DISTANCESpassnp:$jobid_distances_passnp\n" >> $output_dir/logs/jobids.txt
	
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

# Or local	
else
	for count in `seq 1 $sample_count`
	do
		echo "Running stats on sample $count"
		coverage=$( $coverage_cmd $count)

		echo "Running fastqc on sample $count"
        	bamutil_pre=$( $bamutil_preduplicates $count)
	
		if [ $duplicate_filter == "YES" ]; then
                        bamutil_post=$( $bamutil_postduplicates $count)
                fi

	done
	r_script=$( $r_script_cmd)
	distances_filsnp=$( $distances_filsnp_cmd)
	distances_passnp=$( $distances_passnp_cmd)

fi
