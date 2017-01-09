#' @export
svqsc_annotate_QUAL <- function(inputvcffile, model, outputvcffile) {
	# load VCF
	# save original QUAL scores to INFO field svqsc OQ
		# add new INFO to header
		# write QUAL
	# apply model
	# write new QUAL scores to VCF
	# write new VCF
}

#' Generates a calibrated variant scoring model
#'
#' Separate models are used for insertions, deletions, inversions, tandem duplications
#' and inter-chromosomal events as the event signatures for these events can differ
#' significantly.
#'
#' @param vcfgr breakpoints used to train the model
#' @param vcfdf INFO fields reported by the structural variant caller used make the structural variant calls
#' @param modelTransform nrow-invariant function that takes vcfdf and returns
#' only the columns to be used to generate the model. All columns must be numeric.
#' The following additional fields are also available: svLen
#' @param truthgr breakpoints considered true positives.
#' @param requiredSupportingReads truthgr support count required to be called as a true positive.
#' When comparing against a known truth call set, use 1.
#' When comparing directly against long read split alignments, a value greater than 1 is appropriate.
#' @param allowsPartialHits count truth matches at the event level.
#' If TRUE, a variant is called as a true positive if supported by requiredSupportingReads in truthgr
#' If FALSE, a variant is called as a true positive if all constituent breakpoints are supported by requiredSupportingReads in truthgr
#' @return Variant scoring model for use in \code{\link{svqsc_score}}
#' @export
svqsc_train <- function(vcfgr, vcfdf, modelTransform, truthgr, requiredSupportingReads=1, ...) {
	vcfdf <- svqsc_annotate(vcfgr, vcfdf)
	vcfdf <- .svqsc_annotate_tp(vcfgr, vcfdf, truthgr, countOnlyBest, requiredSupportingReads, allowsPartialHits, ...)
	modeldf <- modelTransform(modeldf)
	modeldf$tp <- vcfdf$tp
	modeldf$eventType <- vcfdf$eventType

	model <- .svqsc_generate_model(modeldf, vcfdf, ...)
	model$transform <- modelTransform
	return(model)
}
svqsc_generate_plots <- function(vcfgr, vcfdf, model, truthgr, requiredSupportingReads=1, ...) {
	pred <- svqsc_predict(vcfgr, vcfdf, model)
	vcfdf <- svqsc_annotate(vcfgr, vcfdf)
	vcfdf <- .svqsc_annotate_tp(vcfgr, vcfdf, truthgr, countOnlyBest, requiredSupportingReads, allowsPartialHits, ...)
	modeldf <- modelTransform(modeldf)
	modeldf$tp <- vcfdf$tp
	modeldf$eventType <- vcfdf$eventType
	plots <- list()
	for (event in .eventTypes) {
		eventmodeldf <- modeldf %>%
			dplyr::filter(eventType==event)
			dplyr::select(-eventType)
		plots[[event]] <- list()
		plots[[event]]$pairs <- ggpairs(eventmodeldf %>% mutate(tp=ifelse(tp, "TP", "FP")), columns=2:length(eventmodeldf), mapping=ggplot2::aes(colour=tp, alpha=0.2))
		plots[[event]]$precision_recall <- ggplot() + aes(x=tp, y=precision) +
			# QUAL score is not necessarily in the model so we need to pull it from vcfdf
			geom_line(data=.toPrecRecall(vcfdf$QUAL[vcfdf$eventType==event], eventmodeldf$tp), aes(colour="caller"))
			geom_line(data=.toPrecRecall(pred, eventmodeldf$tp), aes(colour="model"))
	}
	return(plots)
}
.svqsc_annotate_tp <- function(vcfgr, vcfdf, truthgr, requiredSupportingReads, ...) {
	hitscounts <- countBreakpointOverlaps(vcfgr, truthgr, ...)
	if (allowsPartialHits) {
		hitdf <- data.frame(name=names(vcfgr), hitscounts=hitscounts) %>%
			dplyr::group_by(name) %>%
			dplyr::summarise(hitscounts=sum(hitscounts)) %>%
			dplyr::mutate(tp=hitscounts >= requiredSupportingReads) %>%
			dplyr::filter(tp)
	} else {
		hitdf <- data.frame(name=names(vcfgr), hitscounts=hitscounts,
			tp=hitscounts >= requiredSupportingReads) %>%
			dplyr::group_by(name) %>%
			dplyr::summarise(tp=all(tp), hitscounts=sum(hitscounts)) %>%
			dplyr::filter(tp)
	}
	vcfdf$tp <- FALSE
	vcfdf[vcfgr[hitdf$name]$vcfId,]$tp <- TRUE
	return(vcfdf)
}
.svqsc_generate_model <- function(modeldf, vcfdf) {
	assert_that(all(c("eventType", "tp") %in% names(modeldf)))
	# Generate a different model for each event type
	model <- list()
	for (event in .eventTypes) {
		eventmodeldf <- modeldf %>% dplyr::filter(eventType==event)
		if (nrow(eventmodeldf) > 0) {
			cv <- cv.glmnet(eventmodeldf %>% dplyr::select(-tp) %>% as.matrix(), eventmodeldf$tp, alpha=1, family='binomial')
			# TODO: calibrate model probabilities
			model[[event]] <- cv
		}
	}
	return(model)
}
svqsc_annotate <- function(vcfgr, vcfdf) {
	assert_that(all(row.names(vcfdf) %in% vcfgr$vcfId))
	vcfdf$svLen <- 0
	vcfdf[vcfgr$vcfId,]$svLen <- abs(vcfgr$svLen) + abs(vcfgr$insLen)
	vcfdf$eventType <- NA_character_
	vcfdf[vcfgr$vcfId,]$eventType <- .eventType(vcfgr)
	vcfdf$QUAL <- NA_real_
	vcfdf[vcfgr$vcfId,]$QUAL <- vcfgr$QUAL
}
.eventTypes <- c("INS", "DEL", "DUP","INV", "XCHR")
.eventType <- function(vcfgr) {
	p <- partner(vcfgr)
	type <- ifelse(seqnames(vcfgr) != seqnames(p), "XCHR",
		ifelse(strand(vcfgr) == strand(p), "INV",
			ifelse((strand(vcfgr) == "+" & start(vcfgr) > start(p)) | (strand(p) == "+" & start(p) > start(vcfgr)), "DUP",
				ifelse(vcfgr$insLen > vcfgr$svLen / 2, "INS", "DEL"))))
	return(type)
}
#' Score the given variants according to the trained model
#' @param model scoring model generated from \code{\link{svqsc_train}}
svqsc_predict <- function(vcfgr, vcfdf, model) {
	modeldf <- model$transform(vcfdf)
	modeldf <- svqsc_annotate(vcfgr, modeldf)
	modeldf$prediction <- rep(NA_real_, length(vcfgr))
	for (event in .eventTypes) {
		if (!is.null(model[[event]])) {
			eventmodeldf <- modeldf[modeldf$eventType==event,]
			pred <- predict(model[[event]], newx=eventmodeldf %>% as.matrix(), type="response", s="lambda.1se")[,1]
			modeldf[modeldf$eventType==event]$prediction <- pred
		}
	}
	return(modeldf$prediction)
}
.toPrecRecall <- function(scores, tps, rocSlicePoints=NULL) {
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

