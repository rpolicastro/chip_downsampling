#!/usr/bin/env Rscript

library("getopt")
library("Rsubread")

## command line arguments

options <- matrix(c(
	"workdir", "d", 1, "character", "working directory",
	"bam", "b", 1, "character", "original bam file",
	"paired", "e", 1, "character", "if paired end or not",
	"from", "f", 1, "character", "starting value",
	"to", "t", 1, "character", "end value",
	"by", "y", 1, "character", "by value"
), byrow=TRUE, ncol=5)

opt <- getopt(options)
print(opt)

###########################
## Annotating BAM Fragments
###########################

## getting sampling info

samples <- seq(as.numeric(from), as.numeric(to), as.numeric(by))

## preparing peaks file as annotation

options(scipen=999)

peaks <- list()
for (sample in samples) {
	peak.file <- list.files(file.path(workdir, "results", "sampled_peaks"), pattern=paste0(sample,"_.+narrowPeak"))
	peak.data <- read.delim(file.path(workdir, "results", "sampled_peaks", peak.file), header=FALSE, sep="\t")
	peak.data <- peak.data[,c(4,1,2,3,6)]
	colnames(peak.data) <- c("GeneID", "Chr", "Start", "End", "Strand")
	peak.data$Strand <- "*"
	peaks[[as.character(sample)]] <- peak.data
}

## annotating bam fragments to ChEC-seq peaks

dir.create(file.path(opt$workdir, "results", "bam_annotation"))

fractions <- data.frame()
for (sample in samples) {
	bam.file <- file.path(workdir, "results", "sampled_bams", paste0(sample, "_sampled_", basename(opt$bam))) 
	
	# annotating fragments
	if (opt$paired == "TRUE") {
		annotation.results <- featureCounts(
			bam.file,
			annot.ext=peaks[[as.character(sample)]],
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
			bam.file,
			annot.ext=peaks[[as.character(sample)]],
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

	# getting annotation stats
	annotation.stats <- annotation.results$stat
	colnames(annotation.stats) <- c("stats", "counts")
	annotation.stats$frac <- annotation.stats$count / sum(annotation.stats$count)
	
	# adding stats to master data.frame
	df.export <- data.frame(
		"sampling"=sample,
		"peaks"=nrow(peaks[[as.character(sample)]]),
		"overlapping.fragments.count"=annotation.stats[1,2],
		"overlapping.fragments.frac"=annotation.stats[1,3]
	)
	fractions <- rbind(fractions, df.export)
	
	# going back and exporting the raw annotation stats
	outname <- basename(opt$bam)
	outname <- substr(outname, 1, nchar(outname)-4)
	outname <- paste0(sample, "_sampled_", outname, ".tsv")

	write.table(
		annotation.stats, file.path(opt$workdir, "results", "bam_annotation", outname),
		sep="\t", col.names=TRUE, row.names=FALSE, na="", quote=FALSE
	)
}

# exporting data from master data.frame
write.table(
	fractions, file.path(workdir, "results", "bam_annotation", "sampling_results.tsv"),
	sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na=""
)
