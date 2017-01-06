source("human-wgs-50x-common.R", chdir=TRUE)

caller <- "delly"
vcf <- readVcf(paste0("C:/dev/svqsc.data/", caller, ".vcf.gz"), "hg19")
vcfgr <- breakpointRanges(vcf)
seqlevelsStyle(vcfgr) <- "UCSC"
vcfgr <- filterCalls(vcfgr, ignore.altContigs=ignore.altContigs, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, blacklistgr=wgEncodeDacMapabilityConsensusExcludable)
vcfgr$eventType <- eventType(vcfgr)
vcfdf <- unpack(vcf)
vcfdf<- vcfdf[row.names(vcfdf) %in% vcfgr$vcfId,]
#vcfdf <- addEncodeDacColumn(vcfdf, vcfgr)

vcfgr$hits <- countBreakpointOverlaps(vcfgr, truthgr, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest)
vcfgr$tp <- vcfgr$hits >= requiredSupportingReads
vcfdf$tp <- FALSE
vcfdf[vcfgr[vcfgr$tp]$vcfId,]$tp <- TRUE
vcfdf$svLen <- 0
vcfdf[vcfgr$vcfId,]$svLen <- abs(vcfgr$svLen) + abs(vcfgr$insLen)
vcfdf$eventType <- "tmpplaceholder"
vcfdf[vcfgr$vcfId,]$eventType <- vcfgr$eventType
vcfdf$QUAL <- NA_integer_
vcfdf[vcfgr$vcfId,]$QUAL <- vcfgr$QUAL
vcfdf$QUAL <- ifelse(!is.na(vcfdf$QUAL), vcfdf$QUAL, vcfdf$PE + vcfdf$SR)

assert_that(all(row.names(vcfdf) %in% vcfgr$vcfId))

plotPairs <- FALSE
for (event in c("INS", "DEL", "DUP","INV")) {
	subsetdf <- vcfdf %>%
		dplyr::filter(eventType==event)
	if (nrow(subsetdf) > 0) {
		modeldf <- subsetdf %>% dplyr::mutate(logsvlen=log(svLen+1)) %>%
			dplyr::mutate(logPE=log(PE+1), logSR=log(SR+1)) %>%
			dplyr::select(tp, logsvlen, logPE, MAPQ, logSR, SRQ)

		if (plotPairs) {
			# pair-wise comparison
			plot <- ggpairs(modeldf %>% mutate(tp=ifelse(tp, "TP", "FP")), columns=2:length(modeldf), mapping=ggplot2::aes(colour=tp, alpha=0.2))
			ggsave(plot=plot, filename=paste0("plots/", caller, "-pairs-", event, ".png"), units="cm", height=29.7-2, width=21-2)
		}
		# create model
		cv <- cv.glmnet(modeldf %>% dplyr::select(-tp) %>% as.matrix(), modeldf$tp, alpha=1, family='binomial')
		pred <- predict(cv, newx=modeldf %>% dplyr::select(-tp) %>% as.matrix(), type="response", s="lambda.1se")
		# compare
		plot <- ggplot() + aes(x=tp, y=precision) +
			geom_line(data=toPrecRecall(subsetdf$QUAL, subsetdf$tp), aes(colour="caller")) +
			geom_line(data=toPrecRecall(pred[,1], subsetdf$tp), aes(colour="model"))
		ggsave(plot=plot, filename=paste0("plots/", caller, "-roc-", event, ".png"), units="cm", height=29.7-2, width=21-2)
	}
}
library(caret) #install.packages('caret', dependencies=TRUE)
x <- as.matrix(trainingdf %>% select(-tp))
y <- trainingdf$tp
trc = trainControl(method="cv", number=10)
fitM = train(x, y, trControl=trc, method="glmnet", family="binomial", metric="Accuracy")



