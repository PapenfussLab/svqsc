source("human-wgs-50x-functions.R", chdir=TRUE)

datadir <- "C:/dev/svqsc.data/"

# parameters used
mineventsize <- 51 # use the dbSNP/dbVar recommended threshold for minimum SV size
maxeventsize <- 1000000000
ignore.altContigs <- TRUE # ignore alternate contigs
ignore.interchromosomal <- TRUE # we are dealing with germline events
maxgap <- 200 # add error margin for imprecise callers such as breakdancer
ignore.strand <- TRUE # breakdancer does not report direction so this needs to be ignored
sizemargin <- 0.25 # allow +-25% event size
countOnlyBest <- TRUE # count hits only against the best matching variant call (ie consider duplicate variant calls as false)
allowsPartialHits <- FALSE

# plotting
outputdir <- "plots"
theme_set(theme_bw())
suppressWarnings(dir.create(outputdir))


# use the ENCODE DAC blacklist to filter problematic regions
if (!exists("blacklistgr")) {
	wgEncodeDacMapabilityConsensusExcludable <- memoized.load.bed("C:/dev/sv_benchmark/input.na12878/wgEncodeDacMapabilityConsensusExcludable.bed")
	blacklistgr <- wgEncodeDacMapabilityConsensusExcludable
}
if (!exists("truthset")) {
	truthset <- list()
	truthset[["NA12878"]] <- list(
		description="1000 Genomes PacBio", #"1000 Genomes Moleculo", #
		requiredSupportingReads=5, # require at least 5 long reads
		gr=memoized.load.filtered.bedpe(
			paste0(datadir, "/NA12878/truth/", "NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam.sv.bedpe.gz"), # "NA12878.moleculo.bwa-mem.20140110.bam.sv.bedpe.gz"),
			ignore.altContigs, ignore.interchromosomal, mineventsize * (1-sizemargin), maxeventsize * (1+sizemargin), blacklistgr)
	)
	truthset[["HG002"]] <- list(
		description="Genome in a Bottle PacBio",
		requiredSupportingReads=5, # require at least 5 long reads
		gr=memoized.load.filtered.bedpe(
			paste0(datadir, "/HG002/truth/", "HG002_NA24385_son_PacBio_MtSinai_NIST_CSHL_bwamem_bam_GRCh37.bedpe.gz"),
			ignore.altContigs, ignore.interchromosomal, mineventsize * (1-sizemargin), maxeventsize * (1+sizemargin), blacklistgr)
	)
	# Temp hack to adjust bounds until BEDPE is regenerated
	start(truthset[["NA12878"]]$gr) <- start(truthset[["NA12878"]]$gr) - 1
	start(truthset[["HG002"]]$gr) <- start(truthset[["HG002"]]$gr) - 1
}


