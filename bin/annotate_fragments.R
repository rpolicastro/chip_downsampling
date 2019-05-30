#!/usr/bin/env Rscript

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

peaks <- read.delim(opt$peaks, header=TRUE, sep="\t")
peaks <- peaks[,c(4,1,2,3,6)]
colnames(peaks) <- c("GeneID", "Chr", "Start", "End", "Strand")
peaks$Strand <- "*"

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

annotation.stats <- annotation.results$stat
colnames(annotation.stats) <- c("stats", "counts")
annotation.stats$frac <- annotation.stats$count / sum(annotation.stats$count)

outname <- basename(opt$bam)
outname <- substr(outname, 1, nchar(outname)-4)
outname <- paste0(outname, ".tsv")

write.table(
	annotation.stats, file.path(opt$workdir, "results", "bam_annotation", outname),
	sep="\t", col.names=TRUE, row.names=FALSE, na="", quote=FALSE
) 
