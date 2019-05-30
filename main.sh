#!/bin/bash

source settings.conf

#####################
## RIF Like Analysis
#####################

## Preparing Run
## -------------

## activating environment

source active chip-downsampling

## loading settings

source settings.conf

## Downsampling BAM
## ----------------

##  get total number of reads in BAM

READS=$(samtools flagstat $BAM | awk 'NR==1' | cut -d" " -f1)

## divide number by 2 if paired

[[ $PAIRED = TRUE ]] && READS=$(python -c "print $READS/2")

## get sampling percentage

SAMPLE_FRAC=$(python -c "print float($SAMPLE)/float($READS)")

## sample BAM

mkdir -p ${WORKDIR}/results/sampled_bams

samtools view -bs 0$SAMPLE_FRAC $BAM | \
samtools sort -O BAM -o ${WORKDIR}/results/sampled_bams/${SAMPLE}_sampled_$(basename $BAM)

samtools index ${WORKDIR}/results/sampled_bams/${SAMPLE}_sampled_$(basename $BAM)

## Peak Calling
## ------------

mkdir -p ${WORKDIR}/results/sampled_peaks

if [ $PAIRED = TRUE ]
then
	macs2 callpeak \
	-t ${WORKDIR}/results/sampled_bams/${SAMPLE}_sampled_$(basename $BAM) \
	-c $CONTROL \
	-f BAMPE \
	-g $GENOME_SIZE \
	--outdir ${WORKDIR}/results/sampled_peaks \
	-n ${SAMPLE}_sampled_$(basename $BAM .bam)
else
	macs2 callpeak \
	-t ${WORKDIR}/results/sampled_bams/${SAMPLE}_sampled_$(basename $BAM) \
	-c $CONTROL \
	-f BAM \
	-g $GENOME_SIZE \
	--outdir ${WORKDIR}/results/sampled_peaks \
	-n ${SAMPLE}_sampled_$(basename $BAM .bam)
fi

## Annotating BAM fragments against Peaks
## --------------------------------------

Rscript ./bin/annotate_fragments.R \
-d ${WORKDIR} \
-p ${WORKDIR}/results/sampled_peaks/${SAMPLE}_sampled_$(basename $BAM .bam)_peaks.narrowPeak \
-b ${WORKDIR}/results/sampled_bams/${SAMPLE}_sampled_$(basename $BAM) \
-e $PAIRED


