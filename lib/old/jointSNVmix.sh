#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

mkdir $RESULTS_DIR/jointSNVmix/

cd $RESULTS_DIR/jointSNVmix/

date
echo -e "Starting Variant Calling workflow with JointSNVMix\n"

for f in `ls $RESULTS_DIR/Alignment/ | grep "^.*$Control-wodup.bam$"` 
do
	echo -e "4) Variant Analysis"

	filename=`basename $f`
	sample=$(perl $installDir/requisites/processII.pl $filename $Control)

	mkdir $filename
	cd $filename

	bamcontrol="$sample$Control-wodup.bam"
	echo $bamcontrol
	bamcase="$sample$Case-wodup.bam"
	echo $bamcase

    python $JOINT_PATH/jsm.py train $min_normal_depth $skip_size $model $min_tumour_depth $REF_PATH $RESULTS_DIR/Alignment/$bamcontrol $RESULTS_DIR/Alignment/$bamcase parameter_$filename.cfg

    python $JOINT_PATH/jsm.py classify --parameters_file=parameter_$filename.cfg --out_file=$filename.joint --post_process $somatic_threshold $REF_PATH $RESULTS_DIR/Alignment/$bamcontrol $RESULTS_DIR/Alignment/$bamcase 

    #Meter estos valors de 0.8 en el USERConfig
    awk '(($10+$11) > 0.8) && ($18 > 0.8)' $filename.joint > $filename.Somatic.joint
    awk '(($12+$14) > 0.8)' $filename.joint > $filename.LOH.joint
    perl $installDir/requisites/JointToVcf.pl $REF_PATH $filename.LOH.joint
    perl $installDir/requisites/JointToVcf.pl $REF_PATH $filename.Somatic.joint

	vcftools --vcf $filename.Somatic.joint.vcf --out $filename-Joint_ExomeFilter --bed $EXOME_ENRICHMENT --recode --recode-INFO-all
	
	$installDir/requisites/variant_effect_predictor/variant_effect_predictor.pl -i $filename.LOH.joint.vcf --format vcf --output_file $filename.LOH.Joint_Effect.vcf --everything -cache -dir $cacheDir --offline --vcf --force_overwrite
	$installDir/requisites/parserExcell.pl -d $filename.LOH.Joint_Effect.vcf -s jointsnvmix

	$installDir/requisites/variant_effect_predictor/variant_effect_predictor.pl -i $filename.Somatic.joint.vcf --format vcf --output_file $filename.Somatic.Joint_Effect.vcf --everything -cache -dir $cacheDir --offline --vcf --force_overwrite  
	$installDir/requisites/parserExcell.pl -d $filename.Somatic.Joint_Effect.vcf -s jointsnvmix

	cd ..

done
