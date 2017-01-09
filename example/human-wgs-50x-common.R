library(devtools)
#install_github("d-cameron/StructuralVariantAnnotation")
#install_github("d-cameron/svqsc")
library(dplyr)
library(ggplot2)
library(GGally)
library(StructuralVariantAnnotation)
library(assertthat)
library(R.cache)

# parameters used
mineventsize <- 51 # use the dbSNP/dbVar recommended threshold for minimum SV size
maxeventsize <- 1000000000
ignore.altContigs <- TRUE # ignore alternate contigs
ignore.interchromosomal <- TRUE # we are dealing with germline events
maxgap <- 200 # add error margin for imprecise callers such as breakdancer
ignore.strand <- TRUE # breakdancer does not report direction so this needs to be ignored
sizemargin <- 0.25 # allow +-25% event size
countOnlyBest <- TRUE # consider duplicate variant calls as false
considerDuplicateCallsTrue <- FALSE
allowsPartialHits <- FALSE

theme_set(theme_bw())

# use the ENCODE DAC blacklist to filter problematic regions
if (!exists("blacklistgr")) {
	wgEncodeDacMapabilityConsensusExcludable <- memoized.load.bed("C:/dev/sv_benchmark/input.na12878/wgEncodeDacMapabilityConsensusExcludable.bed")
	blacklistgr <- wgEncodeDacMapabilityConsensusExcludable
}
if (!exists("truthgr")) {
	requiredSupportingReads <- 5 # require at least 5 long reads
	truthsample <- "na12878.moleculo"
	#truthgr <- bedpe2breakpointgr("C:/dev/sv_benchmark/input.na12878/lumpy-Mills2012-call-set.bedpe")
	#truthgr <- load.filtered.bedpe("C:/dev/input.na12878/longread/pacbio/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.noid.sv.bedpe.gz",
	truthgr <- memoized.load.filtered.bedpe("C:/dev/input.na12878/longread/moleculo/NA12878.moleculo.clean.bedpe.gz",
		ignore.altContigs,
		ignore.interchromosomal,
		mineventsize * (1-sizemargin),
		maxeventsize * (1+sizemargin),
		blacklistgr)
}


