
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
	result$vcf <- load.filtered.vcf("C:/dev/svqsc.data/", sample, "/", caller, ".vcf.gz", ignore.altContigs, ignore.interchromosomal, mineventsize, maxeventsize, blacklistgr)
	result$vcfgr <- breakpointRanges(result$vcf)
	seqlevelsStyle(result$vcfgr) <- "UCSC"
	result$vcfgr <- filterCalls(result$vcfgr, ignore.altContigs, ignore.interchromosomal, mineventsize, maxeventsize, blacklistgr)
	result$vcfdf <- unpack(result$vcf)
	result$vcfdf <- result$vcfdf[row.names(result$vcfdf) %in% result$vcfgr$vcfId,]
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
	vcf <- vcf[unique(gr$EVENT),]
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
memoized.load.filtered.bedpe <- addMemoization(load.bed)
memoized.load.bed <- addMemoization(load.bed)
memoized.load.svcalls <- addMemoization(load.svcalls)

generate.plots <- function(outputdir="plots", testdfgr, model, truthgr, requiredSupportingReads, maxgap, ignore.strand, sizemargin, countOnlyBest, considerDuplicateCallsTrue, allowsPartialHits) {
	plots <- svqsc_generate_plots(testdfgr$vcfgr, testdfgr$vcfdf, model, truthgr[[train$sample]], requiredSupportingReads, maxgap, ignore.strand, sizemargin, countOnlyBest, considerDuplicateCallsTrue, allowsPartialHits)
	suppressWarnings(dir.create(outputdir))
	for (event in names(model)) {
		prefix <- paste0(outputdir, "/", testdfgr$sample, "-", testdfgr$model_sample,  "-", train$caller, "-", train$event, "-")
		ggsave(plot=plots[[event]]$pairs, filename=paste0(prefix, event, "-pairs.png"), units="cm", height=29.7-2, width=21-2)
		ggsave(plot=plots[[event]]$precision_recall, paste0(prefix, event, "-precision_recall.png"), units="cm", height=29.7-2, width=21-2)
	}
}

