source("human-wgs-50x-common.R", chdir=TRUE)

caller <- "delly"
sample <- "NA12878"
vcf <- readVcf(paste0("C:/dev/svqsc.data/", sample, "/", caller, ".vcf.gz"), "hg19")
vcfgr <- breakpointRanges(vcf)
seqlevelsStyle(vcfgr) <- "UCSC"
vcfgr <- filterCalls(vcfgr, ignore.altContigs=ignore.altContigs, ignore.interchromosomal=ignore.interchromosomal, mineventsize=mineventsize, maxeventsize=maxeventsize, blacklistgr=wgEncodeDacMapabilityConsensusExcludable)
vcfdf <- unpack(vcf)
vcfdf <- vcfdf[row.names(vcfdf) %in% vcfgr$vcfId,]

# Delly doesn't report QUAL scores - default to ranking based on read count
vcfgr$QUAL <- (vcfdf[vcfgr$vcfId,] %>% mutate(score=PE + SR))$score

modelTransform <- function(modeldf) modeldf %>%
	dplyr::mutate(logsvlen=log(svLen+1)) %>%
	dplyr::mutate(logPE=log(PE+1), logSR=log(SR+1)) %>%
	dplyr::select(tp, logsvlen, logPE, MAPQ, logSR, SRQ)

model <- svqsc_train(vcfgr, vcfdf, modelTransform, requiredSupportingReads=requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, considerDuplicateCallsTrue=considerDuplicateCallsTrue, allowsPartialHits=allowsPartialHits)
model_sample <- sample

saveRDS(model, paste("model", caller, "human-wgs-50x", model_sample, ".Rdata", sep="-"))

plots <- svqsc_generate_plots(vcfgr, vcfdf, model, truthgr, requiredSupportingReads=requiredSupportingReads, maxgap=maxgap, ignore.strand=ignore.strand, sizemargin=sizemargin, countOnlyBest=countOnlyBest, considerDuplicateCallsTrue=considerDuplicateCallsTrue, allowsPartialHits=allowsPartialHits)
suppressWarnings(dir.create("plots"))
for (event in names(z)) {
	prefix <- paste0("plots/", sample, "-", model_sample,  "-", caller, "-", event, "-")
	ggsave(plot=plots[[event]]$pairs, filename=paste0(prefix, event, "-pairs.png"), units="cm", height=29.7-2, width=21-2)
	ggsave(plot=plots[[event]]$precision_recall, paste0(prefix, event, "-precision_recall.png"), units="cm", height=29.7-2, width=21-2)
}

# TODO: generate plots for HG002 calls
