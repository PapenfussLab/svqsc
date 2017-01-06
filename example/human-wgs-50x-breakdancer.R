source("human-wgs-50x-common.R", chdir=TRUE)

caller <- "breakdancer"
vcf <- readVcf(paste0("C:/dev/svqsc.data/", caller, ".vcf.gz"), "hg19")
vcfgr <- breakpointRanges(vcf)
seqlevelsStyle(vcfgr) <- "UCSC"
vcfgr <- filterCalls(vcfgr, ignore.altContigs=ignore.altContigs, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, blacklistgr=wgEncodeDacMapabilityConsensusExcludable)
vcfdf <- unpack(vcf)
vcfdf<- vcfdf[row.names(vcfdf) %in% vcfgr$vcfId,]
assert_that(length(vcfgr) == nrow(vcfdf))
#vcfdf <- addEncodeDacColumn(vcfdf, vcfgr)

vcfdf$hits <- countBreakpointOverlaps(vcfgr, truthgr, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest)
vcfdf$tp <- vcfdf$hits >= requiredSupportingReads


modeldf <- vcfdf %>%
	# Orient1 & Orient2 could be used for a strand bias metric
	dplyr::select(IMPRECISE, SVLEN, SVTYPE, Size, Score, num_Reads)


# look at the effect of each column independently
ggpairs(vcfdf. )



result <- list()

breakdancervcf <- readVcf("C:/dev/svqsc/breakdancer.vcf.gz", "hg19")
info(breakdancervcf)[!is.na(info(breakdancervcf)$EVENT)]$Chr1

result[["breakdancer"]] <- go(breakdancervcf, dftransform=function(df, gr) df %>%
		dplyr::select(IMPRECISE, SVLEN, Size, Score, num_Reads))

dellyvcf <- readVcf("C:/dev/svqsc/delly.vcf.gz", "hg19")
fixed(dellyvcf)$QUAL <- (info(dellyvcf)$PE %na% 0) + (info(dellyvcf)$SR %na% 0)
result[["delly"]] <- go(dellyvcf, dftransform=function(df, gr) df %>%
		dplyr::select(PE, MAPQ, SR, SRQ, IMPRECISE, CIPOS))

socratesvcf <- readVcf("C:/dev/svqsc/socrates.vcf.gz", "hg19")
result[["socrates"]] <- go(socratesvcf, dftransform=function(df, gr) df %>%
		# todo mutate(ANCHCONSLEN=nchar(ANCHCONS), REALNCONSLEN=nchar(REALNCONS))
		dplyr::select(QUAL, NLSC, NSSC, BLSC, BSSC, LSSC))

# debug
remove(vcf, df, gr)
considerDuplicateCallsTrue <- FALSE
requiredSupportingReads <- 3
allowsPartialHits <- TRUE
intrachromosomalOnly <- TRUE
traininggr <- gr
trainingdf <- df
hitscounts <- .svqsc_long_read_hits(traininggr, truthgr, considerDuplicateCallsTrue, maxgap=maxgap, ignore.strand=ignore.strand)
trainingdf <- .svqsc_annotate_tp(traininggr, trainingdf, truthgr, considerDuplicateCallsTrue, requiredSupportingReads, allowsPartialHits, maxgap=maxgap, ignore.strand=ignore.strand)
vcf <- socratesvcf

cv <- cv.glmnet(as.matrix(trainingdf %>% select(-tp)), trainingdf$tp, alpha=1, family='binomial')
pred <- predict(cv, newx=as.matrix(df), type="response", s="lambda.1se")

library(caret) #install.packages('caret', dependencies=TRUE)
x <- as.matrix(trainingdf %>% select(-tp))
y <- trainingdf$tp
trc = trainControl(method="cv", number=10)
fitM = train(x, y, trControl=trc, method="glmnet", family="binomial", metric="Accuracy")



