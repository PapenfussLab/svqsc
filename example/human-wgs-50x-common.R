library(devtools)
#install_github("d-cameron/StructuralVariantAnnotation")
#install_github("d-cameron/svqsc")
library(dplyr)
library(ggplot2)
library(GGally)
library(StructuralVariantAnnotation)
library(assertthat)
library(pROC)

theme_set(theme_bw())

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
eventOverlaps <- function(bedgr, vcfdf, vcfgr, ...) {
	vcfdf$overlaps <- FALSE
	vcfdf[vcfgr$vcfId[overlapsAny(vcfgr, bedgr, ...)],]$overlaps <- TRUE
	return(vcfdf$overlaps)
}
addEncodeDacColumn <- function(vcfdf, vcfgr, ...) {
	vcfdf$EncodeDac <- eventOverlaps(wgEncodeDacMapabilityConsensusExcludable, vcfdf, vcfgr, ...)
	return(vcfdf)
}
eventType <- function(vcfgr) {
	p <- partner(vcfgr)
	type <- ifelse(seqnames(vcfgr) != seqnames(p), "ITX",
		ifelse(strand(vcfgr) == strand(p), "INV",
			ifelse((strand(vcfgr) == "+" & start(vcfgr) > start(p)) | (strand(p) == "+" & start(p) > start(vcfgr)), "DUP",
				ifelse(vcfgr$insLen > vcfgr$svLen / 2, "INS", "DEL"))))
	return(type)
}
toPrecRecall <- function(scores, tps, rocSlicePoints=NULL) {
	rocdf <- data.frame(QUAL=scores, tp=tps) %>%
		dplyr::mutate(fp=!tp) %>%
		dplyr::group_by(QUAL) %>%
		dplyr::summarise(tp=sum(tp), fp=sum(fp)) %>%
		dplyr::arrange(dplyr::desc(QUAL)) %>%
		dplyr::mutate(tp=cumsum(tp), fp=cumsum(fp)) %>%
		dplyr::mutate(precision=tp/(tp+fp), fdr=1-precision)
	# calculate area under precision-recall curve
	#rocdf <- rocdf %>%
	#	dplyr::mutate(recall=tp/max(tp))

	if (!is.null(rocSlicePoints)) {
		# subsample along tp and tp+fp axis
		rocdf <- rocdf %>%
			dplyr::slice(unique(c(
				1,
				findInterval(seq(0, max(tp), max(tp)/rocSlicePoints), tp),
				findInterval(seq(0, max(tp + fp), max(tp + fp)/rocSlicePoints), tp + fp),
				n()
			)))
	}
	return(rocdf)
}

# parameters used
mineventsize <- 51 # use the dbSNP/dbVar recommended threshold for minimum SV size
maxeventsize <- 1000000000
ignore.altContigs <- TRUE # ignore alternate contigs
ignore.interchromosomal <- TRUE # we are dealing with germline events
maxgap <- 200 # add error margin for imprecise callers such as breakdancer
ignore.strand <- TRUE # breakdancer does not report direction so this needs to be ignored
sizemargin <- 0.25 # allow +-25% event size
countOnlyBest <- TRUE # consider duplicate variant calls as false
if (!exists("wgEncodeDacMapabilityConsensusExcludable")) {
	# use the ENCODE DAC blacklist to filter problematic regions
	wgEncodeDacMapabilityConsensusExcludable <- import("C:/dev/sv_benchmark/input.na12878/wgEncodeDacMapabilityConsensusExcludable.bed")
	seqlevelsStyle(wgEncodeDacMapabilityConsensusExcludable) <- "UCSC"
}
#truthgr <- bedpe2breakpointgr("C:/dev/sv_benchmark/input.na12878/lumpy-Mills2012-call-set.bedpe")
if (!exists("truthgr")) {
	requiredSupportingReads <- 5 # require at least 5 long reads
	truthgr <- bedpe2breakpointgr("C:/dev/input.na12878/longread/moleculo/NA12878.moleculo.clean.bedpe.gz")
	#truthgr <- bedpe2breakpointgr("C:/dev/input.na12878/longread/pacbio/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.noid.sv.bedpe.gz")
	seqlevelsStyle(truthgr) <- "UCSC"
	tmp <- truthgr
	truthgr <- filterCalls(truthgr,
		ignore.altContigs,
		ignore.interchromosomal,
		mineventsize * (1-sizemargin),
		maxeventsize * (1+sizemargin),
		wgEncodeDacMapabilityConsensusExcludable)
}
