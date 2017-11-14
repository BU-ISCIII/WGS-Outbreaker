#!/bin/bash
## Author S. Monzon
## version v2.0

# Help
# usage: run_variantCalling.sh ....

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

## Folder creation
echo -e "Creating $output_dir/variant_calling"
mkdir -p $output_dir/variant_calling

if [ "$use_sge" = "1" -a $duplicate_filter == "YES" ]; then
 	jobid=$(cat $output_dir/logs/jobids.txt | grep -w "PICARD" | cut -d ':' -f2 )
 	calling_args="${SGE_ARGS} -pe openmp $threads -l h_vmem=$vmem -hold_jid $jobid"
else
 	jobid=$(cat $output_dir/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )
 	calling_args="${SGE_ARGS} -pe openmp $threads -l h_vmem=$vmem -hold_jid $jobid"
fi


if [ $duplicate_filter == "YES" ]; then
	input_list=$duplicateBamArray_list
else
	input_list=$mappingArray_sorted_list
fi

joint_vcf_cmd="$SCRIPTS_DIR/gatk_jointvcf.sh \
        $threads \
        $output_dir/variant_calling/variants_gatk \
        $samples \
        $vcfArray_list \
        $vcffilArray_list \
        $vcfsnpsArray_list \
        $vcfsnpsfilArray_list \
        $vcfindelsArray_list \
        $vcfindelsfilArray_list \
        $ref_path \
        $gatk_path \
        $haplotypeGVCF_list\
        $vcfsnpPass
	$max_snp\
	$window_size"



if [ $know_snps == "NO" ];then
	
	calling_cmd="$SCRIPTS_DIR/gatk_haploid.sh \
	$threads \
	$output_dir/Alignment/BAM  \
	$output_dir/variant_calling/variants_gatk \
	$samples \
	$input_list \
	$vcfArray_list \
	$ref_path \
	$gatk_path \
	$haplotypeGVCF_list"
	

	if [ $variant_calling == "YES" ]; then
                if [ "$use_sge" = "1" ]; then
                calling=$( qsub $calling_args -t 1-$sample_count -N $JOBNAME.CALLING $calling_cmd)
                jobid_calling=$( echo $calling | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "Variant_Calling:$jobid_calling \n" >> $output_dir/logs/jobids.txt
        	
		jointvcf_args="${SGE_ARGS} -pe openmp $threads -l h_vmem=$vmem -hold_jid $jobid_calling"
		jointvcf=$( qsub $jointvcf_args -N $JOBNAME.JOINT.VCF $joint_vcf_cmd)
                jobid_jointvcf=$( echo $jointvcf | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "Joint_ VCF:$jobid_jointvcf \n" >> $output_dir/logs/jobids.txt        
	
		else
		for count in `seq 1 $sample_count`
                do
                        echo "Running variant calling on sample $count"
                        calling=$($calling_cmd $count)
			jointvcf=$(joint_vcf_cmd $count)
                done
                fi
	fi


else
	precalling_cmd="$SCRIPTS_DIR/gatk_preprocessing.sh \
        $threads \
        $output_dir/Alignment/BAM \
        $output_dir/variant_calling/variants_gatk \
        $samples \
        $know_snps \
        $know_indels \
        $input_list \
        $realignedBamArray_list \
        $recalibratedBamArray_list \
        $ref_path \
        $gatk_path"
	
  	calling_cmd="$SCRIPTS_DIR/gatk_haploid.sh \
	$threads \
	$output_dir/variant_calling/variants_gatk/recalibration \
	$output_dir/variant_calling/variants_gatk \
	$samples \
	$recalibratedBamArray_list \
	$vcfsnpsArray_list \
        $vcfsnpsfilArray_list \
        $vcfindelsArray_list \
        $vcfindelsfilArray_list \
        $vcffilArray_list \
        $ref_path \
        $gatk_path \
	$haplotypeGVCF_list \
	$vcfsnpPass"

	if [ $variant_calling == "YES" ]; then
                if [ "$use_sge" = "1" ]; then
                precalling=$( qsub $precalling_args -t 1-$sample_count -N $JOBNAME.CALLING $precalling_cmd)
                jobid_precalling=$( echo $precalling | cut -d ' ' -f3 | cut -d '.' -f1 )
                
		calling_args="${SGE_ARGS} -pe openmp $threads -l h_vmem=$vmem -hold_jid $jobid_precalling"
                calling=$( qsub $calling_args -t 1-$sample_count -N $JOBNAME.CALLING $calling_cmd)
                jobid_calling=$( echo $calling | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "Variant Calling:$jobid_precalling - $jobid_calling \n" >> $output_dir/logs/jobids.txt
        
		jointvcf_arg=calling_args="${SGE_ARGS} -pe openmp $threads -l h_vmem=$vmem -hold_jid $jobid_calling"
                jointvcf=$( qsub $calling_args -N $JOBNAME.JOINT.VCF $joint_vcf_cmd)
                jobid_jointvcf=$( echo $jointvcf | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "Joint vcf:$jobid_jointvcf \n" >> $output_dir/logs/jobids.txt
        
		else
                for count in `seq 1 $sample_count`
                do
                        echo "Running variant calling on sample $count"
                        precalling=$($precalling_cmd $count)
                done
                calling=$($calling_cmd)
      		fi
	fi
fi

vcf_to_tsv_filsnp_cmd="$SCRIPTS_DIR/vcf_to_tsv.sh \
	$output_dir/variant_calling/variants_gatk/variants \
	$vcfsnpsfilArray_list \
	$tsv_filsnp_file"

vcf_to_tsv_passnp_cmd="$SCRIPTS_DIR/vcf_to_tsv.sh \
        $output_dir/variant_calling/variants_gatk/variants \
        $vcfsnpPass \
        $tsv_passnp_file"

tsv_to_msa_filsnp_cmd="$SCRIPTS_DIR/tsv_to_msa.sh \
        $output_dir/variant_calling/variants_gatk/variants \
	$tsv_filsnp_file
	$msa_filsnp_file"

tsv_to_msa_passnp_cmd="$SCRIPTS_DIR/tsv_to_msa.sh \
        $output_dir/variant_calling/variants_gatk/variants \
        $tsv_passnp_file
        $msa_passnp_file"

if [ $vcf_to_msa == "YES" ]; then
	if [ "$use_sge" = "1" ]; then
	vcf_to_tsv_filsnp_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_jointvcf"
	run_vcf_to_tsv_filsnp=$( qsub $vcf_to_tsv_filsnp_args -N $JOBNAME.VCF_TO_TSVfilsnp $vcf_to_tsv_filsnp_cmd)
        jobid_vcf_to_tsv_filsnp=$( echo $run_vcf_to_tsv_filsnp | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "VCF_TO_TSVfilsnp:$jobid_vcf_to_tsv_filsnp\n" >> $output_dir/logs/jobids.txt

	tsv_to_msa_filsnp_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_vcf_to_tsv_filsnp"
        run_tsv_to_msa_filsnp=$( qsub $tsv_to_msa_filsnp_args -N $JOBNAME.TSV_TO_MSAfilsnp $tsv_to_msa_filsnp_cmd)
        jobid_tsv_to_msa_filsnp=$( echo $run_tsv_to_msa_filsnp | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "TSV_TO_MSAfilsnp:$jobid_tsv_to_msa_filsnp\n" >> $output_dir/logs/jobids.txt

	vcf_to_tsv_passnp_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_jointvcf"
        run_vcf_to_tsv_passnp=$( qsub $vcf_to_tsv_passnp_args -N $JOBNAME.VCF_TO_TSVpassnp $vcf_to_tsv_passnp_cmd)
        jobid_vcf_to_tsv_passnp=$( echo $run_vcf_to_tsv_passnp | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "VCF_TO_TSVpassnp:$jobid_vcf_to_tsv_passnp\n" >> $output_dir/logs/jobids.txt

        tsv_to_msa_passnp_args="${SGE_ARGS} -pe openmp $threads -hold_jid $jobid_vcf_to_tsv_passnp"
        run_tsv_to_msa_passnp=$( qsub $tsv_to_msa_passnp_args -N $JOBNAME.TSV_TO_MSApassnp $tsv_to_msa_passnp_cmd)
        jobid_tsv_to_msa_passnp=$( echo $run_tsv_to_msa_passnp | cut -d ' ' -f3 | cut -d '.' -f1 )
        echo -e "TSV_TO_MSApassnp:$jobid_tsv_to_msa_passnp\n" >> $output_dir/logs/jobids.txt


	else

	echo "Running variant calling on sample"
	run_vcf_to_tsv_filsnp=$($vcf_to_tsv_filsnp_cmd)
	run_tsv_to_msa_filsnp=$($tsv_to_msa_filsnp_cmd)
	
	run_vcf_to_tsv_passnp=$($vcf_to_tsv_passnp_cmd)
        run_tsv_to_msa_passnp=$($tsv_to_msa_passnp_cmd)

	fi
fi
