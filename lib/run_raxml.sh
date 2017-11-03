#!/bin/bash
## Author: A. Hernandez
# Help 
# usage: run_raxml.sh ....
#

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

CONFIG_FILE=$1

#Folder creation

mkdir -p $output_dir/RAXML
echo "Directory for RAXML created"

jobid_cfsan_matrix_distance=$(cat $output_dir/logs/jobids.txt | grep -w "MATRIX_DISTANCE_FILES" | cut -d ':' -f2 )
jobid_tsv_to_msa_allsnp=$(cat $output_dir/logs/jobids.txt | grep -w "TSV_TO_MSAallsnp" | cut -d ':' -f2 )
jobid_tsv_to_msa_fil=$(cat $output_dir/logs/jobids.txt | grep -w "TSV_TO_MSAallsnp" | cut -d ':' -f2 )

if [ $cfsan == "YES" ]; then
	mkdir -p $output_dir/RAXML/CFSAN/preser
	mkdir -p $output_dir/RAXML/CFSAN/all_snp
	dir_cfsan=$output_dir/CFSAN
	
	run_raxmlinfe_cfsan_allsnp_cmd="$SCRIPTS_DIR/raxml_inference.sh \
                $dir_cfsan \
                $output_dir/RAXML/CFSAN/all_snp \
                $cfsan_snpma_fasta \
                $model_raxml"

        run_raxmlboot_cfsan_allsnp_cmd="$SCRIPTS_DIR/raxml_bootstrap.sh \
                $dir_cfsan \
                $output_dir/RAXML/CFSAN/all_snp \
                $cfsan_snpma_fasta \
                $model_raxml"

        run_raxmlannot_cfsan_allsnp_cmd="$SCRIPTS_DIR/raxml_annot.sh \
                $output_dir/RAXML/CFSAN/all_snp \
                $model_raxml"


        run_raxmlinfe_cfsan_preser_cmd="$SCRIPTS_DIR/raxml_inference.sh \
                $dir_cfsan \
                $output_dir/RAXML/CFSAN/preser \
                $cfsan_snpma_fil_fasta \
                $model_raxml"

        run_raxmlboot_cfsan_preser_cmd="$SCRIPTS_DIR/raxml_bootstrap.sh \
                $dir_cfsan \
                $output_dir/RAXML/CFSAN/preser \
                $cfsan_snpma_fil_fasta \
                $model_raxml"

        run_raxmlannot_cfsan_preser_cmd="$SCRIPTS_DIR/raxml_annot.sh \
                $output_dir/RAXML/CFSAN/preser \
                $model_raxml"


	if [ "$use_sge" = "1" ]; then
		raxml_cfsan_args="${SGE_ARGS} -hold_jid ${jobid_cfsan_matrix_distance} -pe orte 100 mpirun"
	
#CFSAN all snp

	        raxmlinfe_cfsan_allsnp_qsub=$( qsub -N $JOBNAME.RAXMLinfe_CFSANallsnp $raxml_cfsan_args $run_raxmlinfe_cfsan_allsnp_cmd)
		 jobid_raxmlinfe_cfsan_allsnp=$(echo $raxmlinfe_cfsan_allsnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
       		 echo -e "RAXMLinfe_CFSANallsnp:$jobid_raxmlinfe_cfsan_allsnp\n" >> $output_dir/logs/jobids.txt


        	raxmlboot_cfsanallsnp_args="${SGE_ARGS} -hold_jid ${jobid_raxmlinfe_cfsan_allsnp} -pe orte 100 mpirun"
       		raxmlboot_cfsan_allsnp_qsub=$( qsub -N $JOBNAME.RAXMLboot_CFSANallsnp $raxmlboot_cfsanallsnp_args $run_raxmlboot_cfsan_allsnp_cmd)
		jobid_raxmlboot_cfsan_allsnp=$(echo $raxmlboot_cfsan_allsnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
       		 echo -e "RAXMLboot_CFSANallsnp:$jobid_raxmlboot_cfsan_allsnp\n" >> $output_dir/logs/jobids.txt


        	raxmlannot_cfsan_args="${SGE_ARGS} -hold_jid ${jobid_raxmlboot_cfsan_allsnp}"
        	raxmlannot_cfsan_allsnp_qsub=$( qsub -N $JOBNAME.RAXMLannot_CFSANallsnp $raxmlannot_cfsan_args $run_raxmlannot_cfsan_allsnp_cmd)
       		jobid_raxmlannot_cfsan_allsnp=$(echo $raxmlannot_cfsan_allsnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
       		echo -e "RAXMLannot_CFSANallsnp:$jobid_raxmlannot_cfsan_allsnp\n" >> $output_dir/logs/jobids.txt

#CFSAN preser snp

        	raxmlinfe_cfsan_preser_qsub=$( qsub -N $JOBNAME.RAXMLinfe_CFSANpreser $raxml_cfsan_args $run_raxmlinfe_cfsan_preser_cmd)
        	jobid_raxmlinfe_cfsan_preser=$(echo $raxmlinfe_cfsan_preser_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "RAXMLinfe_CFSANpreser:$jobid_raxmlinfe_cfsan_preser\n" >> $output_dir/logs/jobids.txt


        	raxmlboot_cfsanpreser_args="${SGE_ARGS} -hold_jid ${jobid_raxmlinfe_cfsan_preser} -pe orte 100 mpirun"
        	raxmlboot_cfsan_preser_qsub=$( qsub -N $JOBNAME.RAXMLboot_CFSANpreser $raxmlboot_cfsanpreser_args $run_raxmlboot_cfsan_preser_cmd)
       		jobid_raxmlboot_cfsan_preser=$(echo $raxmlboot_cfsan_preser_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
       		echo -e "RAXMLboot_CFSANpreser:$jobid_raxmlboot_cfsan_preser\n" >> $output_dir/logs/jobids.txt


        	raxmlannot_cfsanpreser_args="${SGE_ARGS} -hold_jid ${jobid_raxmlboot_cfsan_preser}"
        	raxmlannot_cfsan_preser_qsub=$( qsub -N $JOBNAME.RAXMLannot_CFSANpreser $raxmlannot_cfsanpreser_args $run_raxmlannot_cfsan_preser_cmd)
        	jobid_raxmlannot_cfsan_preser=$(echo $raxmlannot_cfsan_preser_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "RAXMLannot_CFSANpreser:$jobid_raxmlannot_cfsan_preser\n" >> $output_dir/logs/jobids.txt

        else

		run_raxmlinfe_cfsan_allsnp=$($run_raxmlinfe_cfsan_allsnp_cmd)
		run_raxmlboot_cfsan_allsnp=$($run_raxmlboot_cfsan_allsnp_cmd)
		run_raxmlannot_cfsan_allsnp=$($run_raxmlannot_cfsan_allsnp_cmd)

		run_raxmlinfe_cfsan_preser=$($run_raxmlinfe_cfsan_preser_cmd)
		run_raxmlboot_cfsan_preser=$($run_raxmlboot_cfsan_preser_cmd)
		run_raxmlannot_cfsan_preser=$($run_raxmlannot_cfsan_preser_cmd)        	


	fi


fi

if [ $variant_calling == "YES" ]; then
	mkdir -p $output_dir/RAXML/GATK/preser
        mkdir -p $output_dir/RAXML/GATK/all_snp
        dir_gatk=$output_dir/variant_calling/variants_gatk/variants/

	run_raxmlinfe_gatk_allsnp_cmd="$SCRIPTS_DIR/raxml_inference.sh \
                $dir_gatk \
                $output_dir/RAXML/GATK/all_snp \
                $msa_allsnp_file \
                $model_raxml"

        run_raxmlboot_gatk_allsnp_cmd="$SCRIPTS_DIR/raxml_bootstrap.sh \
                $dir_gatk \
                $output_dir/RAXML/GATK/all_snp \
                $msa_allsnp_file \
                $model_raxml"

        run_raxmlannot_gatk_allsnp_cmd="$SCRIPTS_DIR/raxml_annot.sh \
                $output_dir/RAXML/GATK/all_snp \
                $model_raxml"

        run_raxmlinfe_gatk_preser_cmd="$SCRIPTS_DIR/raxml_inference.sh \
                $dir_gatk \
                $output_dir/RAXML/GATK/preser \
                $msa_fil_file \
                $model_raxml"

        run_raxmlboot_gatk_preser_cmd="$SCRIPTS_DIR/raxml_bootstrap.sh \
                $dir_gatk \
                $output_dir/RAXML/GATK/preser \
                $msa_fil_file \
                $model_raxml"

        run_raxmlannot_gatk_preser_cmd="$SCRIPTS_DIR/raxml_annot.sh \
                $output_dir/RAXML/GATK/preser \
                $model_raxml"	

	if [ "$use_sge" = "1" ]; then
	        raxml_gatk_args_allsnp="${SGE_ARGS} -hold_jid ${jobid_tsv_to_msa_allsnp} -pe orte 100 mpirun"
		raxml_gatk_args_fil="${SGE_ARGS} -hold_jid ${jobid_tsv_to_msa_fil} -pe orte 100 mpirun"
	
#GATK all snp
        	raxmlinfe_gatk_allsnp_qsub=$( qsub -N $JOBNAME.RAXMLinfe_GATKallsnp $raxml_gatk_args_allsnp $run_raxmlinfe_gatk_allsnp_cmd)
        	jobid_raxmlinfe_gatk_allsnp=$(echo $raxmlinfe_gatk_allsnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "RAXMLinfe_GATKallsnp:$jobid_raxmlinfe_gatk_allsnp\n" >> $output_dir/logs/jobids.txt

        	raxmlboot_gatk_all_args="${SGE_ARGS} -hold_jid ${jobid_raxmlinfe_gatk_allsnp} -pe orte 100 mpirun"
		raxmlboot_gatk_allsnp_qsub=$( qsub -N $JOBNAME.RAXMLboot_GATKallsnp $raxmlboot_gatk_all_args $run_raxmlboot_gatk_allsnp_cmd)
        	jobid_raxmlboot_gatk_allsnp=$(echo $raxmlboot_gatk_allsnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "RAXMLboot_GATKallsnp:$jobid_raxmlboot_gatk_allsnp\n" >> $output_dir/logs/jobids.txt

		raxmlannot_gatk_allsnp_arg="${SGE_ARGS} -hold_jid ${jobid_raxmlboot_gatk_allsnp}"
        	raxmlannot_gatk_allsnp_qsub=$( qsub -N $JOBNAME.RAXMLannot_GATKallsnp $raxmlannot_gatk_allsnp_arg $run_raxmlannot_gatk_allsnp_cmd)
        	jobid_raxmlannot_gatk_allsnp=$(echo $raxmlannot_gatk_allsnp_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "RAXMLannot_GATKallsnp:$jobid_raxmlannot_gatk_allsnp\n" >> $output_dir/logs/jobids.txt

#GATK snp fil

        	raxmlinfe_gatk_preser_qsub=$( qsub -N $JOBNAME.RAXMLinfe_GATKpreser $raxml_gatk_args_fil $run_raxmlinfe_gatk_preser_cmd)
        	jobid_raxmlinfe_gatk_preser=$(echo $raxmlinfe_gatk_preser_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "RAXMLinfe_GATKpreser:$jobid_raxmlinfe_gatk_preser\n" >> $output_dir/logs/jobids.txt

        	raxmlboot_gatk_preser_arg="${SGE_ARGS} -hold_jid ${jobid_raxmlinfe_gatk_preser} -pe orte 100 mpirun"
		raxmlboot_gatk_preser_qsub=$( qsub -N $JOBNAME.RAXMLboot_GATKpreser $raxmlboot_gatk_preser_arg $run_raxmlboot_gatk_preser_cmd)
        	jobid_raxmlboot_gatk_preser=$(echo $raxmlboot_gatk_preser_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "RAXMLboot_GATKpreser:$jobid_raxmlboot_gatk_preser\n" >> $output_dir/logs/jobids.txt

        	raxmlannot_gatk_preser_arg="${SGE_ARGS} -hold_jid ${jobid_raxmlboot_gatk_preser}"
		raxmlannot_gatk_preser_qsub=$( qsub -N $JOBNAME.RAXMLannot_GATKpreser $raxmlannot_gatk_preser_arg $run_raxmlannot_gatk_preser_cmd)
        	jobid_raxmlannot_gatk_preser=$(echo $raxmlannot_gatk_preser_qsub | cut -d ' ' -f3 | cut -d '.' -f1 )
        	echo -e "RAXMLannot_GATKpreser:$jobid_raxmlannot_gatk_preser\n" >> $output_dir/logs/jobids.txt

        else

		run_raxmlinfe_gatk_allsnp=$($run_raxmlinfe_gatk_allsnp_cmd)
                run_raxmlboot_gatk_allsnp=$($run_raxmlboot_gatk_allsnp_cmd)
                run_raxmlannot_gatk_allsnp=$($run_raxmlannot_gatk_allsnp_cmd)

                run_raxmlinfe_gatk_preser=$($run_raxmlinfe_gatk_preser_cmd)
                run_raxmlboot_gatk_preser=$($run_raxmlboot_gatk_preser_cmd)
                run_raxmlannot_gatk_preser=$($run_raxmlannot_gatk_preser_cmd)
	fi

fi
