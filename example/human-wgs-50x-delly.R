source("human-wgs-50x-common.R", chdir=TRUE)

# Use all fields written by DELLY as part of the classifier
modelTransform <- function(modeldf) modeldf %>%
	dplyr::mutate(logsvlen=log(svLen+1)) %>%
	dplyr::mutate(logPE=log(PE+1), logSR=log(SR+1)) %>%
	dplyr::select(logsvlen, logPE, MAPQ, logSR, SRQ)

# Delly doesn't report QUAL scores - default to ranking based on read count
qualtransform <- function(vcfgr, vcfdf) {
	QUAL <- (vcfdf[vcfgr$vcfId,] %>% mutate(score=PE + SR))$score
	return(QUAL)
}

train <- load.svcalls("delly", "NA12878", ignore.altContigs=ignore.altContigs, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, blacklistgr=blacklistgr)
train$vcfdf$QUAL <- qualtransform(train$vcfgr, train$vcfdf)

train$model_sample <- train$sample
train$model <- svqsc_train(train$vcfgr, train$vcfdf, modelTransform, truthgr[[train$sample]], requiredSupportingReads=requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, considerDuplicateCallsTrue=considerDuplicateCallsTrue, allowsPartialHits=allowsPartialHits)

saveRDS(train$model, paste("model", train$caller, "human-wgs-50x", train$model_sample, ".Rdata", sep="-"))

# plot results of training
generate.plots(train, train$model, truthgr[[train$sample]], requiredSupportingReads=requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, considerDuplicateCallsTrue=considerDuplicateCallsTrue, allowsPartialHits=allowsPartialHits)

testdata <- load.svcalls("delly", "HG002", ignore.altContigs=ignore.altContigs, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, blacklistgr=blacklistgr)
testdata$vcfdf$QUAL <- qualtransform(testdata$vcfgr, testdata$vcfdf)

# plot results of training on independent test data
generate.plots(testdata, train$model, truthgr[[testdata$sample]], requiredSupportingReads=requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, considerDuplicateCallsTrue=considerDuplicateCallsTrue, allowsPartialHits=allowsPartialHits)
