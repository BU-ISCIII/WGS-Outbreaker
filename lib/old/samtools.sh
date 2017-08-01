#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

# Creación de las carpetas de resultados del variant calling.
mkdir $RESULTS_DIR/Samtools/
mkdir $RESULTS_DIR/Samtools/Variants
mkdir $RESULTS_DIR/Samtools/ControlVsTumor

# Se van leyendo los ficheros de alineamiento de cada muestra.
for f in `ls $RESULTS_DIR/Alignment/ | grep "wodup.bam$"`
do

filename=`basename $f`
echo -e "3) Variant Calling: $filename"

refname=`basename $REF_PATH`
refdir=`dirname $REF_PATH`

## Se determinan las variantes y el genotipo de cada muestra.
$SCRIPTS_PATH/samtoolsDep pileup -f $REF_PATH -c $RESULTS_DIR/Alignment/$f > $RESULTS_DIR/Samtools/Variants/all_variants_$filename.pileup  & 
$SCRIPTS_PATH/samtoolsDep pileup -f $REF_PATH -cv $RESULTS_DIR/Alignment/$f $minMapQ $minPhredQ -C50 > $RESULTS_DIR/Samtools/Variants/variants_$filename.pileup & 
wait

## Se filtran todas aquellas regiones que no estén en el fichero de enriquecimiento de exomas.
# $SCRIPTS_PATH/pileline/cmd/pileline-rfilter.sh -A $RESULTS_DIR/VariantsSamtools/variants_$filename.pileup ./ -b $EXOME_ENRICHMENT > $RESULTS_DIR/VariantsSamtools/variants_$filename.pileupOntarget 

## Paso de filtrado recomendado por samtools, máxima profundidad 100. Filtrado específico de indels por calidad.
$SCRIPTS_PATH/samtools.pl varFilter -D100 $RESULTS_DIR/Samtools/Variants/variants_$filename.pileup | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' > $RESULTS_DIR/Samtools/Variants/variants_$filename.pileup_filt 

done

## Una vez obtenidas las variantes se realiza la comparación tumor frente a control con PileLine.
 
 for f in `ls $RESULTS_DIR/Samtools/Variants/ | grep "^variants_.*$Control.*.pileup_filt"` 
 do
 	echo -e "4) Variant Analysis"
	
 	## Se procesa el nombre de los ficheros para realizar la comparación.
	filename=`basename $f`
 	sample=$(perl $SCRIPTS_PATH/processII.pl $filename $Control)

 	vcontrol="variants_$sample$Control-wodup.bam.pileup_filt"
 	vcase="variants_$sample$Case-wodup.bam.pileup_filt"

 	allcontrol="all_variants_$sample$Control-wodup.bam.pileup"
 	allcase="all_variants_$sample$Case-wodup.bam.pileup"
 	echo -e "$vcontrol,$vcase,$allcontrol,$allcase"

 	$SCRIPTS_PATH/pileline/cmd/pileline-2smc.sh -a $RESULTS_DIR/Samtools/Variants/$allcontrol -b $RESULTS_DIR/Samtools/Variants/$allcase --variants-a $RESULTS_DIR/Samtools/Variants/$vcontrol --variants-b $RESULTS_DIR/Samtools/Variants/$vcase --out-prefix $RESULTS_DIR/Samtools/ControlVsTumor/$sample 	

 	for f in ` ls $RESULTS_DIR/Samtools/ControlVsTumor/$sample*`
 	do

 	## Finalmente se ejecuta vep para averiguar el posible efecto de las mutaciones obtenidas.
 	filename=`basename $f` 
 	perl $SCRIPTS_PATH/variant_effect_predictor/variant_effect_predictor.pl -i $f --format pileup --output_file $RESULTS_DIR/Samtools/ControlVsTumor/$filename-ef.txt --vcf --everything -cache $cacheDir --offline --force_overwrite

 	done

 done

