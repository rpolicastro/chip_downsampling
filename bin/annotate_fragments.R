#!/usr/bin/env Rscript

library("tidyverse")
library("getopt")
library("Rsubread")

## command line arguments

options <- matrix(c(
	"workdir", "d", 1, "character", "working directory",
	"peaks", "p", 1, "character", "peak file",
	"bam", "b", 1, "character", "bam file",
	"paired", "e", 1, "character", "if paired end or not"
), byrow=TRUE, ncol=5)

opt <- getopt(options)
print(opt)

###########################
## Annotating BAM Fragments
###########################

## preparing peaks file as annotation

peaks <- read.delim(opt$peaks, header=TRUE, sep="\t") %>%
	select(4, 1, 2, 3, 6) %>%
	rename("GeneID"=1, "Chr"=2, "Start"=3, "End"=4, "Strand"=5) %>%
	mutate("Strand"="*")

## annotating bam fragments to ChEC-seq peaks

if (opt$paired == "TRUE") {
	annotation.results <- featureCounts(
		opt$bam,
		annot.ext=peaks,
		isGTFAnnotationFile=FALSE,
		useMetaFeatures=FALSE,
		minOverlap=10,
		largestOverlap=TRUE,
		countMultiMappingReads=FALSE,
		strandSpecific=0,
		isPairedEnd=TRUE,
	)
} else {
        annotation.results <- featureCounts(
                opt$bam,
                annot.ext=peaks,
                isGTFAnnotationFile=FALSE,
                useMetaFeatures=FALSE,
                minOverlap=10,
                largestOverlap=TRUE,
                countMultiMappingReads=FALSE,
                strandSpecific=0,
                isPairedEnd=FALSE,
		readExtension3=200
        )
}

## getting annotation stats

dir.create(file.path(opt$workdir, "results", "bam_annotation"))

annotation.stats <- annotation.results$stat %>%
	rename("stats"=1, "count"=2) %>%
	mutate(frac=count/sum(count))

outname <- basename(opt$bam) %>% substr(., 1, nchar(.)-4) %>% paste0(., ".tsv")
write.table(
	annotation.stats, file.path(opt$workdir, "results", "bam_annotation", outname),
	sep="\t", col.names=TRUE, row.names=FALSE, na="", quote=FALSE
) 
