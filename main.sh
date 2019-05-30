#!/bin/bash

source settings.conf

#####################
## RIF Like Analysis
#####################

## activating environment

source active chip-downsampling

## loading settings

source settings.conf

## downsampling BAM

# get total number of reads in BAM
READS=$(samtools flagstat $BAM | awk 'NR==1' | cut -d" " -f1)

# divide number by 2 if paired
[[ $PAIRED = TRUE ]] && READS=$(python -c "print $READS/2")

# get sampling percentage
SAMPLE_FRAC=$(python -c "print float($SAMPLE)/float($READS)")

# sample BAM
mkdir -p ${WORKDIR}/results

samtools view -bs 0$SAMPLE_FRAC $BAM | \
samtools sort -O BAM -o ${WORKDIR}/results/${SAMPLE}_sampled_$(basename $BAM)

samtools index ${WORKDIR}/results/${SAMPLE}_sampled_$(basename $BAM)
