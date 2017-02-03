#!/usr/bin/env Rscript
library("optparse")

source("human-wgs-50x-common.R", chdir=TRUE)

# Caller Model Transforms
cmt <- list()

# Breakdancer
cmt$breakdancer <- list()
cmt$breakdancer$ReadCount <- function(modeldf) modeldf %>%
  dplyr::mutate(ReadCount=num_Reads) %>%
  dplyr::select(ReadCount)
cmt$breakdancer$Score <- function(modeldf) modeldf %>%
  dplyr::select(Score)
cmt$breakdancer$logSize <- function(modeldf) modeldf %>%
  dplyr::mutate(logSize=log(abs(Size)+1)) %>%
  dplyr::select(logSize)
cmt$breakdancer$logsvlen_lognum_Reads_Size_Score <- function(modeldf) modeldf %>%
	dplyr::mutate(logsvlen=log(svLen+1)) %>%
	dplyr::mutate(lognum_Reads=log(num_Reads+1)) %>%
	dplyr::select(logsvlen, lognum_Reads, Size, Score)

# crest
cmt$crest <- list()
cmt$crest$ReadCount <- function(modeldf) modeldf %>%
  dplyr::mutate(ReadCount=(left_softclipped_read_count %na% 0) + (right_softclipped_read_count %na% 0)) %>%
  dplyr::select(ReadCount)
cmt$crest$allNumeric <-
	dplyr::mutate(logsvlen=log(svLen+1)) %>%
	dplyr::select(logsvlen,
		left_assembly_length %na% 0,   left_softclipped_read_count %na% 0,  left_percent_identity %na% 0,  left_percent_multimapping %na% 0,  left_coverage %na% 0,
		right_assembly_length %na% 0, right_softclipped_read_count %na% 0, right_percent_identity %na% 0, right_percent_multimapping %na% 0, right_coverage %na% 0)

# delly
cmt$delly$ReadCount <- function(modeldf) modeldf %>%
  dplyr::mutate(ReadCount=PE+SR) %>%
  dplyr::select(ReadCount)
cmt$delly$logsvlen_logPE_MAPQ_logSR_SRQ <- function(modeldf) modeldf %>%
	dplyr::mutate(logsvlen=log(svLen+1)) %>%
	dplyr::mutate(logPE=log(PE+1), logSR=log(SR+1)) %>%
	dplyr::select(logsvlen, logPE, MAPQ, logSR, SRQ)

# manta
cmt$manta$QUAL <- function(modeldf) modeldf %>%
  dplyr::select(QUAL)
cmt$manta$numericBreakpoint <- function(modeldf) modeldf %>%
	dplyr::mutate(logsvlen=log(svLen+1)) %>%
  # TODO: what to do about asymetrical breakend fields?
	dplyr::select(logsvlen,
	  QUAL,
	  BND_DEPTH %na% 0,
	  BND_PAIR_COUNT %na% 0,
	  DOWNSTREAM_PAIR_COUNT %na% 0,
	  UPSTREAM_PAIR_COUNT %na% 0,
	  HOMLEN %na% 0,
	  PAIR_COUNT %na% 0,
	  SVINSLEN %na% 0
	  )

# manta
cmt$manta$QUAL <- function(modeldf) modeldf %>%
  dplyr::select(QUAL %na% 0)
cmt$manta$allNumeric <- function(modeldf) modeldf %>%
	dplyr::mutate(logsvlen=log(svLen+1)) %>%
  # TODO: what to do about asymetrical results?
	dplyr::select(logsvlen,
	  BND_DEPTH %na% 0,
	  BND_PAIR_COUNTH %na% 0,
	  DOWNSTREAM_PAIR_COUNT %na% 0,
	  HOMLEN %na% 0,
	  PAIR_COUNT %na% 0,
	  SVINSLEN %na% 0,
	  UPSTREAM_PAIR_COUNT %na% 0)

# gridss
cmt$gridss$QUAL <- function(modeldf) modeldf %>%
  dplyr::select(QUAL)
cmt$gridss$logCountsModel  <- function(modeldf) {
	partnerdf <- modeldf[modeldf$PARID,]
	modeldf <- modeldf %>%
		dplyr::mutate(logsvlen=log(svLen+1)) %>%
		dplyr::mutate(
			logREF=log(REF + partnerdf$REF + 1),
			logREFPAIR=log(REFPAIR + partnerdf$REFPAIR + 1),
			logBSC=log(BSC + partnerdf$BSC + 1),
			logBUM=log(BUM + partnerdf$BUM + 1),
			logRP=log(RP + 1),
			ReciprocalAssembly=AS > 0 & RAS > 0,
			Assembly=AS > 0 | RAS > 0) %>%
		dplyr::mutate(IHOMLEN=abs(IHOMPOS.1)+abs(IHOMPOS.2)) %>%
		dplyr::select(
			logsvlen,
			ReciprocalAssembly, Assembly,
			HOMLEN, IHOMLEN,
			#SR, SRQ, #TODO: why are these zeroes?
			logRP, RPQ,
			QUAL,
			logBSC, logBUM,
			logREF, logREFPAIR)
	return(modeldf)
}

for (caller in names(cmt)) {
  if (!interactive() && !(caller %in% commandArgs(trailingOnly=TRUE))) {
    # check command-line args and only process that data subset
    next
  }
  train <- load.svcalls(caller, "NA12878", ignore.altContigs=ignore.altContigs, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, blacklistgr=blacklistgr)
  testdata <- load.svcalls(caller, "HG002", ignore.altContigs=ignore.altContigs, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, blacklistgr=blacklistgr)

  models <- list()
  for (modelName in names(cmt[[caller]])) {
    models[[modelName]] <- svqsc_train(train$vcfgr, train$vcfdf, cmt[[caller]][[modelName]], truthset[[train$sample]]$gr, requiredSupportingReads=truthset[[train$sample]]$requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, allowsPartialHits=allowsPartialHits)
    models[[modelName]]$name <- paste0(modelName, "_", train$sample)

    saveRDS(models[[modelName]], paste("model", train$caller, "human-wgs-50x", models[[modelName]]$name, ".Rdata", sep="-"))
  }
  # plot results of training
  generate.plots(outputdir=outputdir, testdfgr=train, models, truthgr=truthset[[train$sample]]$gr, requiredSupportingReads=truthset[[train$sample]]$requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, allowsPartialHits=allowsPartialHits)
  # plot results of training on independent test data
  generate.plots(outputdir=outputdir, testdata, models, truthset[[testdata$sample]]$gr, requiredSupportingReads=truthset[[testdata$sample]]$requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, allowsPartialHits=allowsPartialHits)
}

dfpred <- bind_rows(lapply(list.files(path="plots", pattern="*predictions.tsv"), function(filename) {
	read.table(paste0("plots/", filename))
}))
# TODO: add models:
# - raw caller
# - always 0
# - always 1
dfpred %>%
	dplyr::group_by(sample, caller, eventType, modelName) %>%
	dplyr::summarise(
		# https://en.wikipedia.org/wiki/Brier_score
		observed_tp_frequency=sum(tp)/n(),
		brier=sum((prediction-tp)**2)/n())
		# how many probably categories do we want? Deciles?

# calibration or validation curve
ggplot(dfpred %>%
		dplyr::mutate(prediction_category=round(prediction, 1)) %>%
		dplyr::group_by(sample, caller, eventType, modelName, prediction_category) %>%
		dplyr::summarise(prediction=mean(prediction), actual=mean(tp), n=n())) +
	aes(x=prediction, y=actual) +
	geom_point(aes(size=n, colour=sample)) +
	geom_line(data=bind_rows(
		dfpred %>% dplyr::group_by(sample, caller, eventType, modelName) %>% summarise(prediction=0, actual=0),
		dfpred %>% dplyr::group_by(sample, caller, eventType, modelName) %>% summarise(prediction=1, actual=1))) +
	facet_grid(caller + modelName ~ eventType)
# TODO: phred scale of the above





