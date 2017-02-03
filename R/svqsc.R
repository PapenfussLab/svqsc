#@export
#svqsc_annotate_QUAL <- function(inputvcffile, model, outputvcffile) {
	# load VCF
	# save original QUAL scores to INFO field svqsc OQ
		# add new INFO to header
		# write QUAL
	# apply model
	# write new QUAL scores to VCF
	# write new VCF
#}

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
svqsc_train <- function(vcfgr, vcfdf, modelTransform, truthgr,
		countOnlyBest=TRUE,
		requiredSupportingReads=1, allowsPartialHits=FALSE,
		maxgap = 0L, minoverlap = 1L, ignore.strand = FALSE,
		sizemargin = 0.25, restrictMarginToSizeMultiple = 0.5) {
	vcfdf <- svqsc_annotate(vcfgr, vcfdf)
	vcfdf <- svqsc_annotate_truth(vcfgr, vcfdf, truthgr, countOnlyBest=countOnlyBest, requiredSupportingReads=requiredSupportingReads, allowsPartialHits=allowsPartialHits, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, sizemargin=sizemargin, restrictMarginToSizeMultiple=restrictMarginToSizeMultiple)
	modeldf <- modelTransform(vcfdf)
	modeldf$tp <- vcfdf$tp
	modeldf$eventType <- vcfdf$eventType

	model <- .svqsc_generate_model(modeldf, vcfdf)
	model$transform <- modelTransform
	return(model)
}
svqsc_evaluate_model <- function(vcfgr, vcfdf, model, truthgr,
		countOnlyBest=TRUE,
		requiredSupportingReads=1, allowsPartialHits=FALSE,
		maxgap = 0L, minoverlap = 1L, ignore.strand = FALSE,
		sizemargin = 0.25, restrictMarginToSizeMultiple = 0.5) {
  pred <- svqsc_predict(vcfgr, vcfdf, model)
	vcfdf <- svqsc_annotate(vcfgr, vcfdf)
	vcfdf <- svqsc_annotate_truth(vcfgr, vcfdf, truthgr, countOnlyBest=countOnlyBest, requiredSupportingReads=requiredSupportingReads, allowsPartialHits=allowsPartialHits, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, sizemargin=sizemargin, restrictMarginToSizeMultiple=restrictMarginToSizeMultiple)
	modeldf <- model$transform(vcfdf)
	modeldf$eventType <- vcfdf$eventType
	modeldf$tp <- vcfdf$tp
	modeldf$prediction <- pred
	modeldf$QUAL <- vcfdf$QUAL
	result <- list()
	for (event in .eventTypes) {
	  eventmodeldf <- modeldf %>% filter(eventType==event)
		result[[event]] <- list()
		precRecallModel <- .toPrecRecall(pred[vcfdf$eventType==event], eventmodeldf$tp)
		result[[event]]$precision_recall <- precRecallModel
		result[[event]]$precision_recall_auc <- .precRecallAUC(precRecallModel)
	}
	result$model <- model
	result$modeldf <- modeldf
	return(result)
}
svqsc_annotate_truth <- function(vcfgr, vcfdf, truthgr,
		countOnlyBest, requiredSupportingReads, allowsPartialHits,
		maxgap, minoverlap, ignore.strand,
		sizemargin, restrictMarginToSizeMultiple) {
	hitscounts <- countBreakpointOverlaps(vcfgr, truthgr, countOnlyBest=countOnlyBest, maxgap=maxgap, minoverlap=minoverlap, ignore.strand=ignore.strand, sizemargin=sizemargin, restrictMarginToSizeMultiple=restrictMarginToSizeMultiple)
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
#'
#' @param minevents Minimum number of true/false events for a model to be trained
.svqsc_generate_model <- function(modeldf, vcfdf, minevents=16) {
	assert_that(all(c("eventType", "tp") %in% names(modeldf)))
	# Generate a different model for each event type
	model <- list()
	for (event in .eventTypes) {
		eventmodeldf <- modeldf %>%
			dplyr::filter(eventType==event) %>%
			dplyr::select(-eventType)
		# need a minimum number of true and false events to be able to train the model
		if (all(c(sum(eventmodeldf$tp), sum(!eventmodeldf$tp)) >= minevents)) {
			write(paste("Training", event), stderr())
			cv <- cv.glmnet(model.matrix(~., eventmodeldf %>% dplyr::select(-tp)), eventmodeldf$tp, alpha=1, family='binomial')
			# TODO: calibrate model probabilities
			model[[event]] <- cv
		} else {
			write(paste("Not training", event, "due to insufficient data"), stderr())
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
	if (!("QUAL" %in% names(vcfdf))) {
		vcfdf$QUAL <- NA_real_
		vcfdf[vcfgr$vcfId,]$QUAL <- vcfgr$QUAL
	}
	return(vcfdf)
}
.eventTypes <- c("INS", "DEL", "DUP","INV", "XCHR")
.eventType <- function(vcfgr) {
	p <- partner(vcfgr)
	type <- ifelse(seqnames(vcfgr) != seqnames(p), "XCHR",
		ifelse(strand(vcfgr) == strand(p), "INV",
			ifelse((strand(vcfgr) == "+" & start(vcfgr) > start(p)) | (strand(p) == "+" & start(p) > start(vcfgr)), "DUP",
				ifelse(vcfgr$insLen > abs(vcfgr$svLen) / 2, "INS", "DEL"))))
	return(type)
}
#' Score the given variants according to the trained model
#' @param model scoring model generated from \code{\link{svqsc_train}}
svqsc_predict <- function(vcfgr, vcfdf, model) {
	vcfdf <- svqsc_annotate(vcfgr, vcfdf)
	modeldf <- model$transform(vcfdf)
	modeldf$prediction <- rep(NA_real_, nrow(modeldf))
	modeldf$eventType <- vcfdf$eventType
	for (event in .eventTypes) {
		cv <- model[[event]]
		if (!is.null(cv)) {
			eventmodeldf <- modeldf %>%
				dplyr::filter(eventType==event) %>%
				dplyr::select(-eventType, -prediction)
			pred <- predict(cv, newx=model.matrix(~., eventmodeldf), type="response", s="lambda.1se")[,1]
			modeldf[modeldf$eventType==event,]$prediction <- pred
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
.precRecallAUC <- function(rocdf) {
	return(rocdf %>%
		# recall of 1 is the maximum recall of the caller,
		# not a recall of all events
		dplyr::mutate(recall=tp/max(tp)) %>%
		dplyr::mutate(areaToNext=(lead(recall) - recall)*(lead(precision) + precision) / 2) %>%
		dplyr::summarise(auc=sum(ifelse(is.na(areaToNext), 0, areaToNext)))
		)$auc
}

