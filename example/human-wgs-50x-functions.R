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

generate.plots <- function(outputdir, testdfgr, models, truthgr, requiredSupportingReads, maxgap, ignore.strand, sizemargin, countOnlyBest, allowsPartialHits) {
	prefix <- paste0(outputdir, "/", testdfgr$sample, "-", testdfgr$caller, "-")
	evalmodellist <- lapply(models, function(model) {
		result <- svqsc_evaluate_model(testdfgr$vcfgr, testdfgr$vcfdf, model, truthgr, requiredSupportingReads=requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, allowsPartialHits=allowsPartialHits)
	})
	# plot ROC curves
	dfprecrecall <- bind_rows(lapply(evalmodellist, function(em) {
		bind_rows(lapply(.eventTypes, function(et) {
			df <- em[[et]]$precision_recall
			if (!is.null(df) && nrow(df) > 0) {
				df$eventType <- et
				df$modelName <- em$model$name
				return(df)
			}
			return(NULL)
		}))
	}))
	ggplot(dfprecrecall) + aes(x=tp, y=precision, colour=modelName) +
		geom_line() +
		facet_wrap(~ eventType)
		labs(x="True positives", y="Precision", colour="Model")
	ggsave(file=paste0(prefix, "precision_recall.png"), units="cm", height=10, width=21-2)
	# AUC
	dfprecrecallauc <- bind_rows(lapply(evalmodellist, function(em) {
		bind_rows(lapply(.eventTypes, function(et) {
			auc <- em[[et]]$precision_recall_auc
			if (!is.null(auc)) {
				return(data.frame(
					caller=caller,
					eventType=et,
					modelName=em$model$name,
					precision_recall_auc=auc,
				  stringsAsFactors=FALSE))
			}
			return(NULL)
		}))
	}))
	write.table(dfprecrecallauc, file=paste0(prefix, "precision_recall_auc.tsv"))

	# Calculate Brier score
	dfpred <- bind_rows(lapply(evalmodellist, function(em) {
		sample <- testdfgr$sample
		caller <- testdfgr$caller
		return (em$modeldf %>%
			dplyr::mutate(sample=sample, caller=caller, caller_prediction=10**(-QUAL/10), modelName=em$model$name) %>%
			dplyr::select(sample, caller, eventType, modelName, tp, prediction, caller_prediction))
		}))
	write.table(dfpred, file=paste0(prefix, "predictions.tsv"))

	for (result in evalmodellist) {
		model <- result$model
		for (event in names(result)) {
			if (event %in% .eventTypes) {
				tryCatch({
					eventmodeldf <- result$modeldf %>%
					  dplyr::filter(eventType==event) %>%
					  dplyr::select(-eventType, -prediction)
					if (nrow(eventmodeldf) > 0) {
  					pairsplot <- ggpairs(eventmodeldf %>% mutate(tp=ifelse(tp, "TP", "FP")),
  						#columns=1:length(eventmodeldf-2),
  						mapping=ggplot2::aes(colour=tp, alpha=0.2))
  					ggsave(plot=pairsplot, filename=paste0(prefix, event, "-", model$name, "-pairs.png"), units="cm", height=29.7-2, width=21-2)
					}
				})
			}
		}
	}
}

#' Fills in NAs from the values from the partner breakend
fill.na.from <- function(df, partnerId) {
	for (col in names(df)) {
		df[is.na(df[[col]]),][[col]] <- df[df[is.na(df[[col]]),][[partnerId]],][[col]]
	}
	return(df)
}

