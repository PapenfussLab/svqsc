library(devtools)
#install_github("d-cameron/StructuralVariantAnnotation")
#install_github("d-cameron/svqsc")
library(dplyr)
library(ggplot2)
library(GGally)
library(StructuralVariantAnnotation)
library(assertthat)
library(R.cache)

filterCalls <- function(gr, ignore.altContigs, ignore.interchromosomal, mineventsize, maxeventsize, blacklistgr) {
	if (ignore.altContigs) {
		gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y")) & seqnames(partner(gr)) %in% paste0("chr", c(1:22, "X", "Y")),]
	}
	if (is.null(gr$insLen)) {
		# TODO: how do we get insertion length from a bedpe?
		gr$insLen <- 0
	}
	if (is.null(gr$svLen)) {
		gr$svLen <- ifelse(seqnames(gr) == seqnames(partner(gr)), abs(start(gr) - start(partner(gr))) + gr$insLen, NULL)
	}
	if (ignore.interchromosomal) {
		gr <- gr[seqnames(gr) == seqnames(partner(gr))]
	}
	if (!is.null(mineventsize)) {
		gr <- gr[is.na(gr$svLen) | abs(gr$svLen) >= mineventsize]
	}
	if (!is.null(maxeventsize)) {
		gr <- gr[is.na(gr$svLen) | abs(gr$svLen) <= maxeventsize]
	}
	if (!is.null(blacklistgr)) {
		gr <- gr[!overlapsAny(gr, blacklistgr) & !overlapsAny(partner(gr), blacklistgr)]
	}
	assert_that(all(names(partner(gr)) %in% names(gr)))
	return(gr)
}
addEncodeDacColumn <- function(vcfdf, vcfgr, ...) {
	vcfdf$EncodeDac <- eventOverlaps(wgEncodeDacMapabilityConsensusExcludable, vcfdf, vcfgr, ...)
	return(vcfdf)
}
load.svcalls <- function(caller, sample, ignore.altContigs, ignore.interchromosomal, mineventsize, maxeventsize, blacklistgr) {
	result <- list()
	result$caller <- caller
	result$sample <- sample
	filename <- paste0(datadir, "/", sample, "/", caller, ".vcf.gz")
	result$vcf <- load.filtered.vcf(filename, ignore.altContigs, ignore.interchromosomal, mineventsize, maxeventsize, blacklistgr)
	if (nrow(result$vcf) == 0) {
		stop(paste("No VCF records in ", filename))
	}
	result$vcfgr <- breakpointRanges(result$vcf)
	seqlevelsStyle(result$vcfgr) <- "UCSC"
	result$vcfgr <- filterCalls(result$vcfgr, ignore.altContigs, ignore.interchromosomal, mineventsize, maxeventsize, blacklistgr)
	result$vcfdf <- unpack(result$vcf)
	result$vcfdf <- result$vcfdf[row.names(result$vcfdf) %in% result$vcfgr$vcfId,]
	return(result)
}
load.filtered.vcf <- function(file, ignore.altContigs, ignore.interchromosomal, mineventsize, maxeventsize, blacklistgr) {
	vcf <- readVcf(file, "hg19")
	gr <- breakpointRanges(vcf)
	seqlevelsStyle(gr) <- "UCSC"
	gr <- filterCalls(gr,
		ignore.altContigs,
		ignore.interchromosomal,
		mineventsize,
		maxeventsize,
		blacklistgr)
	vcf <- vcf[unique(gr$vcfId),]
	return(vcf)
}
load.filtered.bedpe <- function(file, ignore.altContigs, ignore.interchromosomal, mineventsize, maxeventsize, blacklistgr) {
	gr <- bedpe2breakpointgr(file)
	seqlevelsStyle(gr) <- "UCSC"
	gr <- filterCalls(gr,
		ignore.altContigs,
		ignore.interchromosomal,
		mineventsize,
		maxeventsize,
		blacklistgr)
	return(gr)
}
load.bed <- function(file) {
	bedgr <- import(file)
	seqlevelsStyle(bedgr) <- "UCSC"
	return(bedgr)
}
# R.cache memoization
memoized.load.filtered.vcf <- addMemoization(load.filtered.vcf)
memoized.load.filtered.bedpe <- addMemoization(load.filtered.bedpe)
memoized.load.bed <- addMemoization(load.bed)
memoized.load.svcalls <- addMemoization(load.svcalls)

generate.plots <- function(outputdir, testdfgr, model, truthgr, requiredSupportingReads, maxgap, ignore.strand, sizemargin, countOnlyBest, allowsPartialHits) {
	plots <- svqsc_generate_plots(testdfgr$vcfgr, testdfgr$vcfdf, model, truthgr, requiredSupportingReads=requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, allowsPartialHits=allowsPartialHits)
	for (event in names(model)) {
		if (event %in% .eventTypes) {
			prefix <- paste0(outputdir, "/", testdfgr$sample, "-", testdfgr$caller, "-", event, "-")
			ggsave(plot=plots[[event]]$pairs, filename=paste0(prefix, "pairs.png"), units="cm", height=29.7-2, width=21-2)
			ggsave(plot=plots[[event]]$precision_recall, paste0(prefix, "precision_recall-", model$name,".png"), units="cm", height=29.7-2, width=21-2)
		}
	}
}


