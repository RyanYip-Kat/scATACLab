library(GenomicRanges)
library(Matrix)
library(BSgenome)
library(dplyr)
keepGR<-function(gr,genome){
	if(tolower(genome%in%c("hg38","hg19"))){
		keepChr=paste0("chr",c(1:22,"X"))
	}else{
		keepChr=paste0("chr",c(1:19,"X"))
	}
	gr <- gr[seqnames(gr) %in% keepChr]
        gr <- keepSeqlevels(gr, keepChr)
	gr
}	
callSummitsMACS2 <- function(
	bedFile = NULL,
	pathToMacs2 = "macs2",
	genomeSize = 2.7e9, 
	shift = -75, 
	extsize = 150, 
	cutOff = 0.05, 
	method = "q",
	additionalParams = "--nomodel --nolambda",
	outdir="result",
	prefix="macs2",
	save=TRUE,
	cleanUp=FALSE
	){

	stopifnot(tolower(method) %in% c("p","q"))
	stopifnot(!is.null(genomeSize))
	utility <- .checkPath(pathToMacs2)

	#Output Files
	#bedName <- gsub("\\.insertions.bed", "", bedFile)
	bedName=prefix
	summitsFile <- paste0(bedName, "_summits.bed")
	narrowPeaksFile <- paste0(bedName, "_peaks.narrowPeak")
	xlsFile <- paste0(bedName, "_peaks.xls")

	#Create MACS2 Command
	cmd <- sprintf("callpeak -g %s --name %s --treatment %s --outdir %s --format BED --call-summits --keep-dup all %s", 
		genomeSize, basename(bedName), bedFile, dirname(bedName), additionalParams)

	if(!is.null(shift) & !is.null(extsize)){
		cmd <- sprintf("%s --shift %s --extsize %s", cmd , shift, extsize)
	}

	if(tolower(method) == "p"){
		cmd <- sprintf("%s -p %s", cmd , cutOff)
	}else{
		cmd <- sprintf("%s -q %s", cmd , cutOff)
	}


	run <- system2(pathToMacs2, cmd, wait=TRUE, stdout=NULL, stderr=NULL)

	#Read Summits!
	out <- data.table::fread(summitsFile, select = c(1,2,3,5))
	out <- GRanges(out$V1, IRanges(out$V2 + 1, out$V3), score = out$V5)
	if(cleanUp){
		#Remove Files
	        r2 <- suppressWarnings(file.remove(summitsFile, narrowPeaksFile, xlsFile))
	}
	if(save){
		outFile=file.path(outdir,paste0(prefix,"_summits.rds"))
		saveRDS(out,outFile)
	}
	return(out)

}

.getQuantiles<-function (v = NULL, len = length(v)) 
{
    if (length(v) < len) {
        v2 <- rep(0, len)
        v2[seq_along(v)] <- v
    }
    else {
        v2 <- v
    }
    p <- trunc(rank(v2))/length(v2)
    if (length(v) < len) {
        p <- p[seq_along(v)]
    }
    return(p)
}

nonOverlappingGR<-function (gr = NULL, by = "score", decreasing = TRUE, verbose = FALSE) 
{
    stopifnot(by %in% colnames(mcols(gr)))
    .clusterGRanges <- function(gr = NULL, filter = TRUE, by = "score", 
        decreasing = TRUE) {
        gr <- sort(sortSeqlevels(gr))
        r <- GenomicRanges::reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE)
        o <- findOverlaps(gr, r, ignore.strand = TRUE)
        mcols(gr)$cluster <- subjectHits(o)
        gr <- gr[order(mcols(gr)[, by], decreasing = decreasing), 
            ]
        gr <- gr[!duplicated(mcols(gr)$cluster), ]
        gr <- sort(sortSeqlevels(gr))
        mcols(gr)$cluster <- NULL
        return(gr)
    }
    if (verbose) {
        message("Converging", appendLF = FALSE)
    }
    i <- 0
    grConverge <- gr
    while (length(grConverge) > 0) {
        if (verbose) {
            message(".", appendLF = FALSE)
        }
        i <- i + 1
        grSelect <- .clusterGRanges(gr = grConverge, filter = TRUE, 
            by = by, decreasing = decreasing)
        grConverge <- subsetByOverlaps(grConverge, grSelect, 
            invert = TRUE, ignore.strand = TRUE)
        if (i == 1) {
            grAll <- grSelect
        }
        else {
            grAll <- c(grAll, grSelect)
        }
    }
    message(sprintf("Converged after %s iterations!", i))
    if (verbose) {
        message("\nSelected ", length(grAll), " from ", length(gr))
    }
    grAll <- sort(sortSeqlevels(grAll))
    return(grAll)
}

identifyReproduciblePeaks <- function(
	summitFiles = NULL,
	summitNames = NULL,
	reproducibility = 0.51,
        extendSummits = 250,
        blacklist = NULL,
        prefix = NULL
       ){

  errorList <- mget(names(formals()),sys.frame(sys.nframe()))

	nonOverlapPassES <- tryCatch({

	  summits <- lapply(seq_along(summitFiles), function(x){
	  	grx <- readRDS(summitFiles[x])
	  	grx <- subsetByOverlaps(grx, blacklist, invert = TRUE) #Not Overlapping Blacklist!
	  	grx$GroupReplicate <- paste0(summitNames[x])
	  	grx
	  })
	  summits <- Reduce("c", as(summits, "GRangesList"))

	  extendedSummits <- resize(summits, extendSummits * 2 + 1, "center")
	  extendedSummits <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
	    nonES <- nonOverlappingGR(x, by = "score", decreasing = TRUE)
	    nonES$replicateScoreQuantile <- round(.getQuantiles(nonES$score),3)
	    nonES
	  })
	  extendedSummits <- Reduce("c", as(extendedSummits, "GRangesList"))

	  nonOverlapES <- nonOverlappingGR(extendedSummits, by = "replicateScoreQuantile", decreasing = TRUE)

	  overlapMat <- lapply(split(extendedSummits, extendedSummits$GroupReplicate), function(x){
	    overlapsAny(nonOverlapES, x)
	  }) %>% Reduce("cbind", .)

	  if(length(summitFiles) > 1){
	    nonOverlapES$Reproducibility <- rowSums(overlapMat)
	    nonOverlapES$ReproducibilityPercent <- round(rowSums(overlapMat) / ncol(overlapMat) , 3)
	    n <- length(summitFiles)
	    minRep <- eval(parse(text=reproducibility))
	    if(!is.numeric(minRep)){
	    	stop("Error reproducibility not numeric when evaluated!")
	    }
	  	idxPass <- which(nonOverlapES$Reproducibility >= minRep)
	  	nonOverlapPassES <- nonOverlapES[idxPass]
	  }else{
	    nonOverlapES$Reproducibility <- rep(NA, length(nonOverlapES))
	    nonOverlapPassES <- nonOverlapES
	  }

	  nonOverlapPassES$groupScoreQuantile <- round(.getQuantiles(nonOverlapPassES$replicateScoreQuantile),3)
	  mcols(nonOverlapPassES) <- mcols(nonOverlapPassES)[,c("score","replicateScoreQuantile", "groupScoreQuantile", "Reproducibility", "GroupReplicate")]

	  nonOverlapPassES

	}, error = function(e){
	       message("INFO : identifyReproduciblePeaks Error!")

	})

  return(nonOverlapPassES)

}

.validGRanges<-function (gr = NULL)
{
    stopifnot(!is.null(gr))
    if (inherits(gr, "GRanges")) {
        return(gr)
    }
    else {
        stop("Error cannot validate genomic range!")
    }
}

.checkPath<-function (u = NULL, path = NULL, throwError = TRUE) 
{
    if (is.null(u)) {
        out <- TRUE
    }
    out <- lapply(u, function(x, error = TRUE) {
        if (Sys.which(x) == "") {
            if (!is.null(path) && file.exists(file.path(path, 
                x))) {
                o <- TRUE
            }
            else {
                if (throwError) {
                  stop(x, " not found in path, please add ", 
                    x, " to path!")
                }
                else {
                  o <- FALSE
                }
            }
        }
        else {
            o <- TRUE
        }
        return(o)
    }) %>% unlist %>% all
    return(out)
}


validBSgenome<-function (genome = NULL, masked = FALSE) 
{
    stopifnot(!is.null(genome))
    if (inherits(genome, "BSgenome")) {
        return(genome)
    }
    else if (is.character(genome)) {
        genome <- tryCatch({
            bsg <- eval(parse(text = genome))
            if (inherits(bsg, "BSgenome")) {
                return(bsg)
            }
            else {
                stop("genome is not a BSgenome valid class!")
            }
        }, error = function(x) {
            BSgenome::getBSgenome(genome, masked = masked)
        })
        return(genome)
    }
    else {
        stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
    }
}

extendGR<-function (gr = NULL, upstream = NULL, downstream = NULL)
{
    isMinus <- BiocGenerics::which(strand(gr) == "-")
    isOther <- BiocGenerics::which(strand(gr) != "-")
    start(gr)[isOther] <- start(gr)[isOther] - upstream
    end(gr)[isOther] <- end(gr)[isOther] + downstream
    end(gr)[isMinus] <- end(gr)[isMinus] + upstream
    start(gr)[isMinus] <- start(gr)[isMinus] - downstream
    return(gr)
}

fastAnnoPeaks <- function(
	peaks = NULL,
	BSgenome = NULL,
	geneAnnotation = NULL,
	promoterRegion = c(2000, 100)
	){

	#Validate
	peaks <- .validGRanges(peaks)
	peakSummits <- resize(peaks,1,"center")
	geneAnnotation$genes <- .validGRanges(geneAnnotation$genes)
	geneAnnotation$exons <- .validGRanges(geneAnnotation$exons)
	geneAnnotation$TSS <- .validGRanges(geneAnnotation$TSS)
	BSgenome <- validBSgenome(BSgenome)

	#First Lets Get Distance to Nearest Gene Start
	distPeaks <- distanceToNearest(peakSummits, resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
	mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
	mcols(peaks)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
	promoters <- extendGR(resize(geneAnnotation$genes, 1, "start"), upstream = promoterRegion[1], downstream = promoterRegion[2])
	op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
	og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
	oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
	type <- rep("Distal", length(peaks))
	type[which(og & oe)] <- "Exonic"
	type[which(og & !oe)] <- "Intronic"
	type[which(op)] <- "Promoter"
	mcols(peaks)$peakType <- type

	#First Lets Get Distance to Nearest TSS's
	distTSS <- distanceToNearest(peakSummits, resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
	mcols(peaks)$distToTSS <- mcols(distTSS)$distance
	if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
		mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distPeaks)]
	}else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
		mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distPeaks)]
	}

	#Get NucleoTide Content
	nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks))
  	mcols(peaks)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  	mcols(peaks)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  	peaks

}

get_args<-function(){
	library(argparse)
        parser <- ArgumentParser(description='CallPeaks for BED files ...')
        parser$add_argument("--bedfiles",
		    nargs="+",
                    type="character",
                    default=NULL,
		    help="BED files")

        parser$add_argument("--names",
		    nargs="+",
		    default=NULL,
		    type="character",
		    help="names for BED file according")

        parser$add_argument("--outdir",
                    type="character",
                    default="PeakCall")


        parser$add_argument("--genome",
                    type="character",
                    default="mm10",
                    choices=c("hg38","hg19","mm9","mm10"),
                    help="which genome to be used")
	args <- parser$parse_args()
	return(args)
}
source("/home/ye/Work/R/scATAC/ArchR/utils/CreateGeneAnnotation.R")
args=get_args()
outdir=args$outdir
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=TRUE)
}

############################## paramaters
genome=args$genome
geneAnnotation=createGeneAnnotation(genome)
if(genome=="hg38"){
	library(BSgenome.Hsapiens.UCSC.hg38)
	bsgenome=BSgenome.Hsapiens.UCSC.hg38
	genomeSize=2.7e9
	keepChr=paste0("chr",c(1:22,"X"))
	blacklist=Signac::blacklist_hg38
}else{
	library(BSgenome.Mmusculus.UCSC.mm10)
	bsgenome=BSgenome.Mmusculus.UCSC.mm10
	genomeSize=1.87e9
	keepChr=paste0("chr",c(1:19,"X"))
	blacklist=Signac::blacklist_mm10
}

BEDFiles=args$bedfiles
Names=args$names

############################### callpeaks
message("INFO : CallPeaks ...")
summitFiles<-c()
summitNames<-c()

for(i in seq_along(BEDFiles)){
	name=Names[i]
	bedfile=BEDFiles[i]
	cat(sprintf("INFO :   %d of %d  ---  %s \n",i,length(BEDFiles),name))
	peak=callSummitsMACS2(bedFile=bedfile,
			 genomeSize=genomeSize,
			 outdir=outdir,
			 prefix=name,
			 cleanUp=FALSE,
			 save=TRUE)
	peak <- keepGR(peak,genome)
	#peak <- fastAnnoPeaks(peak,bsgenome,geneAnnotation)
	keepChrPeak=file.path(outdir,paste0(name,"_clean_summits.rds"))
	saveRDS(peak,keepChrPeak)
	summitFiles<-c(summitFiles,keepChrPeak)
	summitNames<-c(summitNames,name)
}

message("INFO : identify ReproduciblePeaks...")
gr=identifyReproduciblePeaks(summitFiles=summitFiles,
			     summitNames=summitNames,
			     blacklist=blacklist,
			     reproducibility = 0.51,
			     extendSummits = 250)
message("INFO : anno Peaks ...")
gr <- fastAnnoPeaks(gr,bsgenome,geneAnnotation)
saveRDS(gr,file.path(outdir,"AnnoPeaks_summits.rds"))
message("INFO : Done!")
