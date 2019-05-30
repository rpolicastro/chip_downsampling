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

READS=$(samtools flagstat $BAM | awk 'NR==1' | cut -d" " -f1)
