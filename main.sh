#!/bin/bash

#####################
## RIF Like Analysis
#####################

## Preparing Run
## -------------

## activating environment

conda activate chip-downsampling

## loading settings

source settings.conf

## Downsampling BAM
## ----------------

##  get total number of reads in BAM

READS=$(samtools flagstat $BAM | awk 'NR==1' | cut -d" " -f1)

## divide number by 2 if paired

[[ $PAIRED = TRUE ]] && READS=$(python -c "print $READS/2")

## get sampling percentages

SAMPLES=($(seq $FROM $BY $TO))

SAMPLE_FRACS=()
for SAMPLE in ${SAMPLES[@]}; do
	SAMPLE_FRACS+=($(python -c "print float($SAMPLE)/float($READS)"))
done

## sample BAM

mkdir -p ${WORKDIR}/results/sampled_bams

N_SAMPLES=${#SAMPLE_FRACS[@]}
for N in $(seq 0 1 $(bc <<< $N_SAMPLES-1)); do
	samtools view -bs 0${SAMPLE_FRACS[$N]} $BAM | \
	samtools sort -O BAM -o ${WORKDIR}/results/sampled_bams/${SAMPLES[$N]}_sampled_$(basename $BAM)

	samtools index ${WORKDIR}/results/sampled_bams/${SAMPLES[$N]}_sampled_$(basename $BAM)
done

## Peak Calling
## ------------

mkdir -p ${WORKDIR}/results/sampled_peaks

for SAMPLE in ${SAMPLES[@]}; do
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
done

## Annotating BAM fragments against Peaks
## --------------------------------------

Rscript ./bin/annotate_fragments.R \
-d ${WORKDIR} \
-p ${WORKDIR}/results/sampled_peaks/${SAMPLE}_sampled_$(basename $BAM .bam)_peaks.narrowPeak \
-b ${WORKDIR}/results/sampled_bams/${SAMPLE}_sampled_$(basename $BAM) \
-e $PAIRED


