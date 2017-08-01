#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

Creación de carpetas de resultados
echo -e "Starting Variant Calling workflow"
date
echo -e "\n"
mkdir $RESULTS_DIR/VariantsGatk
mkdir $RESULTS_DIR/VariantsGatk/Realignment/
mkdir $RESULTS_DIR/VariantsGatk/Recalibration/
mkdir $RESULTS_DIR/VariantsGatk/Variants/
mkdir $RESULTS_DIR/VariantsGatk/OnlyTumor/

for file in $RESULTS_DIR/Alignment/*-wodup.bam    
do
   name=`basename $file`
	echo -e "1) Local Realignment: $file"	

   # Generar información de read en el fichero BAM, es solicitado por GATK.
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/picard/AddOrReplaceReadGroups.jar I=$file O=$file-RG.bam RGID=HWI-ST731_6 RGLB=undefined RGPL=Illumina RGPU=undefined RGSM=$file VALIDATION_STRINGENCY=LENIENT  

   # Indexar el nuevo fichero BAM generado.
   $SAMTOOLS_PATH/samtools index $file-RG.bam

   #java -Xmx2G -jar $SCRIPTS_PATH/picard/ReorderSam.jar I=$file-RG.bam O=$file-RG_sorted.bam REFERENCE=$REF_PATH
   $SAMTOOLS_PATH/samtools sort $file-RG.bam $file-RG_sorted

   $SAMTOOLS_PATH/samtools index $file-RG_sorted.bam	
   
   ## Realineamiento local - preparación.
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
   		-T RealignerTargetCreator \
   		-I $file-RG_sorted.bam \
      -known $INDELS_1000G \
   		-R $REF_PATH \
   		-o $RESULTS_DIR/VariantsGatk/$name-IndelRealigner.intervals \
   		-nt $THREADS \
   		-S LENIENT \
   		-log $RESULTS_DIR/VariantsGatk/$name-targetCreator.log

   ## Realineamiento alrededor de indels.
	java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
   		-T IndelRealigner \
   		-I $file-RG_sorted.bam \
         -known $INDELS_1000G \
   		-R $REF_PATH \
   		-targetIntervals $RESULTS_DIR/VariantsGatk/$name-IndelRealigner.intervals \
   		-o $RESULTS_DIR/VariantsGatk/Realignment/$name-realigned.bam \
   		-S LENIENT \
 #   		-log $RESULTS_DIR/VariantsGatk/$name-realigner.log

   echo -e "2) Base Quality Recalibration"

   ## Recalibración de la calidad de las bases.
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
	   -T BaseRecalibrator \
	   -R $REF_PATH \
   	-I $RESULTS_DIR/VariantsGatk/Realignment/$name-realigned.bam \
      -knownSites $SNPS_hg19 \
      -knownSites $HAPMAP_33 \
   	-cov ReadGroupCovariate \
   	-cov QualityScoreCovariate \
   	-cov CycleCovariate \
   	-cov ContextCovariate \
    -o $RESULTS_DIR/VariantsGatk/Recalibration/$name-recal_data.grp \
   	-nct $THREADS \
   	-S LENIENT \
   	-log $RESULTS_DIR/VariantsGatk/$name-recal.log


    java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R $REF_PATH \
    -I $RESULTS_DIR/VariantsGatk/Realignment/$name-realigned.bam \
    -BQSR $RESULTS_DIR/VariantsGatk/Recalibration/$name-recal_data.grp \
    -o $RESULTS_DIR/VariantsGatk/Recalibration/$name-recalibrated.bam \
    -nct $THREADS \
    -S LENIENT \
    -log $RESULTS_DIR/VariantsGatk/$name-print_recal.log

   echo -e "Variant calling with Unified Genotyper"

   # Variant calling. Sólo variantes.
  java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -T UnifiedGenotyper \
     -R $REF_PATH \
     -I $RESULTS_DIR/VariantsGatk/Recalibration/$name-recalibrated.bam \
      --dbsnp $SNPS_hg19 \
      -o $RESULTS_DIR/VariantsGatk/Variants/$name-snps.raw.vcf \
      -stand_call_conf 30.0 \
      -stand_emit_conf 10.0 \
      -glm BOTH \
      -nct $THREADS \
      -S LENIENT \
      -log $RESULTS_DIR/VariantsGatk/$name-variants.log

   # Sitios variantes e invariantes.
  # java -Djava.io.tmpdir=$TEMP -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
  # -l INFO \
  # -R $REF_PATH \
  # -T UnifiedGenotyper \
  # -I $RESULTS_DIR/VariantsGatk/Recalibration/$name-recalibrated.bam \
  # -L $EXOME_ENRICHMENT \
  # --dbsnp $SNPS_hg19 \
  # -o $RESULTS_DIR/VariantsGatk/Variants/allsites-$name.vcf \
  # --output_mode EMIT_ALL_SITES
  # -nct $THREADS \
  # -S LENIENT \
  # -log $RESULTS_DIR/VariantsGatk/$name-allvariants.log

   # Se seleccionan sólo los snps y se filtran.
   echo -e "Select and Filter SNPs"
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
      -T SelectVariants \
      -R $REF_PATH \
      -V $RESULTS_DIR/VariantsGatk/Variants/$name-snps.raw.vcf \
      -o $RESULTS_DIR/VariantsGatk/Variants/$name-snps.vcf \
      -selectType SNP \
      -nt $THREADS \
      -S LENIENT \
      -log $RESULTS_DIR/VariantsGatk/$name-selectSNPs.log

   # Se filtran según los parámetros por defecto que recomienda GATK. El filtro por cobertura es muy bajo
   # y el verdadero filtro se realizará después de realizar la comparación tumor frente a control. 
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
      -R $REF_PATH \
      -T VariantFiltration \
      -V $RESULTS_DIR/VariantsGatk/Variants/$name-snps.vcf \
      -o $RESULTS_DIR/VariantsGatk/Variants/$name.SNPs.filtered.vcf \
      --clusterWindowSize 10 \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --filterExpression "DP <5 " \
      --filterName "LowCoverage" \
      --filterExpression "QUAL <30.0 " \
      --filterName "VeryLowQual" \
      --filterExpression "QUAL >30.0 && QUAL <50.0 " \
      --filterName "LowQual" \
      --filterExpression "QD <1.5 " \
      --filterName "LowQD" \
      --filterExpression "SB >-10.0 " \
      --filterName "StrandBias" \
      --filterExpression "FS >60.0 " \
      --filterName "p-value StrandBias" \
      --filterExpression "HaplotypeScore >13.0 " \
      --filterName "HaplotypeScore" \
      -S LENIENT \
      -log $RESULTS_DIR/VariantsGatk/$name-filterSNPs.log

   ## Realizamos el mismo procedimiento con los indels.
   echo -e "Select and Filter Indels"
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
      -R $REF_PATH \
      -T SelectVariants \
      -V $RESULTS_DIR/VariantsGatk/Variants/$name-snps.raw.vcf \
      -o $RESULTS_DIR/VariantsGatk/Variants/$name-indels.vcf \
      -selectType INDEL \
      -nt $THREADS \
      -S LENIENT \
      -log $RESULTS_DIR/VariantsGatk/$name-selectIndels.log


   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
      -T VariantFiltration \
      -R $REF_PATH \
      -V $RESULTS_DIR/VariantsGatk/Variants/$name-indels.vcf \
      -o $RESULTS_DIR/VariantsGatk/Variants/$name-indels_filtered.vcf \
      --filterExpression "QD < 2.0 || MQ0 > 50" \
      --filterName "Nov09filters" \
      -S LENIENT \
      -log $RESULTS_DIR/VariantsGatk/$name-filterIndels.log
done

for f in $RESULTS_DIR/Alignment/*$Control-wodup.bam
do
   echo -e "4) Variant Analysis"

   filename=`basename $f`
   sample=$(perl $SCRIPTS_PATH/processII.pl $filename $Control)

   vcontrol="$sample$Control-wodup.bam.SNPs.filtered.vcf"
   vcase="$sample$Case-wodup.bam.SNPs.filtered.vcf"

   ## Se realiza la comparación tumor frente a control utilizando vcftools.

   bgzip $RESULTS_DIR/VariantsGatk/Variants/$vcase ; tabix -p vcf $RESULTS_DIR/VariantsGatk/Variants/$vcase.gz
   bgzip $RESULTS_DIR/VariantsGatk/Variants/$vcontrol; tabix -p vcf $RESULTS_DIR/VariantsGatk/Variants/$vcontrol.gz

   vcf-isec -c -f $RESULTS_DIR/VariantsGatk/Variants/$vcase.gz $RESULTS_DIR/VariantsGatk/Variants/$vcontrol.gz > $RESULTS_DIR/VariantsGatk/OnlyTumor/$sample-Onlytumor.vcf
   
   ## Se filtran las variantes por cobertura y por calidad. Cogiendo sólo aquellas que pasan todos los filtros. AQUÍ SE PUEDE PONER EL FILTRO PARA NOVELS grep \.\
   vcf-annotate -f Q=100/d=10 $RESULTS_DIR/VariantsGatk/OnlyTumor/$sample-Onlytumor.vcf | grep "PASS" > $RESULTS_DIR/VariantsGatk/OnlyTumor/$sample-Onlytumor_FilterPASS.vcf

   ## Se analiza el efecto de las mutaciones encontradas.
   $SCRIPTS_PATH/variant_effect_predictor/variant_effect_predictor.pl -i $RESULTS_DIR/VariantsGatk/OnlyTumor/$sample-Onlytumor_FilterPASS.vcf --format vcf --output_file $RESULTS_DIR/VariantsGatk/OnlyTumor/$filename-GATK_Effect.vcf -cache -dir $cacheDir --everything --offline --vcf --force_overwrite
   perl $SCRIPTS_PATH/parserExcell.pl -d $RESULTS_DIR/VariantsGatk/OnlyTumor/$filename-GATK_Effect.vcf -s gatk
done