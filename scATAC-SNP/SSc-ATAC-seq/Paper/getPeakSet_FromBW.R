library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(rtracklayer)
library(readr)
library(GenomicRanges)
library(magrittr)
library(stringr)
library(edgeR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(argparse)
set.seed(1)

countInsertions <- function(query, fragments,name, by = "score"){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  colnames(sparseM)=paste(name,1:ncol(sparseM),sep="-")
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

extendedPeakSet <- function(df, BSgenome = NULL, extend = 250, blacklist = NULL, nSummits = 100000){
  #Helper Functions
  readSummits <- function(file){
    df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end","name","score"))))
    df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
    return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
  }
  nonOverlappingGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
    stopifnot(by %in% colnames(mcols(gr)))
    clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE){
      gr <- sort(sortSeqlevels(gr))
      r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
      o <- findOverlaps(gr,r)
      mcols(gr)$cluster <- subjectHits(o)
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      gr <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(gr))
      mcols(gr)$cluster <- NULL
      return(gr)
    }
    if(verbose){
      message("Converging", appendLF = FALSE)
    }
    i <-  0
    gr_converge <- gr
    while(length(gr_converge) > 0){
      if(verbose){
        message(".", appendLF = FALSE)
      }
      i <-  i + 1
      gr_selected <- clusterGRanges(gr = gr_converge, filter = TRUE, by = by, decreasing = decreasing)
      gr_converge <- subsetByOverlaps(gr_converge ,gr_selected, invert=TRUE) #blacklist selected gr
      if(i == 1){ #if i=1 then set gr_all to clustered
        gr_all <- gr_selected
      }else{
        gr_all <- c(gr_all, gr_selected)
      } 
    }
    if(verbose){
      message("\nSelected ", length(gr_all), " from ", length(gr))
    }
    gr_all <- sort(sortSeqlevels(gr_all))
    return(gr_all)
  }
  #Check-------
  stopifnot(extend > 0)
  stopifnot("samples" %in% colnames(df))
  stopifnot("groups" %in% colnames(df))
  stopifnot("summits" %in% colnames(df))
  stopifnot(!is.null(BSgenome))
  stopifnot(all(apply(df,1,function(x){file.exists(paste0(x[3]))})))
  #------------
  #Deal with blacklist
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- rtracklayer::import.bed(blacklist)
  }
  stopifnot(inherits(blacklist,"GenomicRanges"))
  #------------
  #Time to do stuff
  chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  groups <- unique(df$groups)
  groupGRList <- GenomicRanges::GenomicRangesList(lapply(seq_along(groups), function(i){
      df_group = df[which(df$groups==groups[i]),]
      grList <- GenomicRanges::GenomicRangesList(lapply(paste0(df_group$summits), function(x){
        extended_summits <- readSummits(x) %>%
          resize(., width = 2 * extend + 1, fix = "center") %>%     
          subsetByOverlaps(.,chromSizes,type="within") %>%
          subsetByOverlaps(.,blacklist,invert=TRUE) %>%
          nonOverlappingGRanges(., by="score", decreasing=TRUE)
        extended_summits <- extended_summits[order(extended_summits$score,decreasing=TRUE)]
        if(!is.null(nSummits)){
          extended_summits <- head(extended_summits, nSummits)
        }
        mcols(extended_summits)$scoreQuantile <- trunc(rank(mcols(extended_summits)$score))/length(mcols(extended_summits)$score)
        extended_summits
      }))
      #Non Overlapping
      grNonOverlapping <- nonOverlappingGRanges(unlist(grList), by = "scoreQuantile", decreasing = TRUE)
      #Free Up Memory
      remove(grList)
      gc()
      grNonOverlapping
    }))
  grFinal <- nonOverlappingGRanges(unlist(groupGRList), by = "scoreQuantile", decreasing = TRUE)
  grFinal <- sort(sortSeqlevels(grFinal))
  return(grFinal)
}

groupSums <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

makedir=function(path){
	if(!dir.exists(path)){
		dir.create(path,recursive=TRUE)
	}
}

########################################   callpeak paramenters
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--meta",
                    type="character",
                    default=NULL,
                    help="*.bw file in meta csv file")


args <- parser$parse_args()


method <- "q"
cutoff <- 0.05
shift <- -75
extsize <- 150
genome_size <- 2.7e9
genome <- BSgenome.Hsapiens.UCSC.hg38
macs2="/home/ye/anaconda3/envs/scatac/bin/macs2"

########################################
message("Write Paper bw file into BED files")
DF=read.csv(args$meta,header=FALSE,sep=",",stringsAsFactors=FALSE)
bwFiles=as.character(DF$V2)
Names=as.character(DF$V1)


outPath=args$outdir
dirPeaks= file.path(outPath,"bw2CallPeaks")
dirBeds=dirPeaks
outDir=file.path(outPath,"Counts")

makedir(outDir)
makedir(dirPeaks)
makedir(dirBeds)
for(i in 1:length(bwFiles)){
   name=Names[i]
   bw=bwFiles[i]
   if(!file.exists(file.path(dirPeaks,paste0(name,"_summits.bed")))){
	#message(paste0("INFO : ", bw))
	cat(sprintf("INFO : [ %d of %d ] --- [ %s ] \n",i,length(bwFiles),bw))
	gr=import(bw)
	bed=as.data.frame(gr)
	out=bed[,c("seqnames","start","end")]

	SampleBED=file.path(dirBeds,paste0(name,".bed"))
	message(paste0("INFO : Write BED File into : ",SampleBED))
	write.table(out,SampleBED,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

	message(sprintf("CallPeak %s from %s",name,SampleBED))
	cmdPeaks <- sprintf(
	    "%s callpeak -g %s --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all",
	    macs2,
	    genome_size,
	    name,
	    SampleBED,
	    dirPeaks
	  )
	if (!is.null(shift) & !is.null(extsize)) {
	  cmdPeaks <- sprintf("%s --shift %s --extsize %s", cmdPeaks, shift, extsize)
	}
	if (tolower(method) == "p") {
	  cmdPeaks <- sprintf("%s -p %s", cmdPeaks, cutoff)
	}else {
	  cmdPeaks <- sprintf("%s -q %s", cmdPeaks, cutoff)
	}
	message("Running Macs2...")
	message(cmdPeaks)
        system(cmdPeaks, intern = TRUE)
	system(paste0("rm ",SampleBED))
   }
}


#-------------------------------------------------------------------------------------------------
# Make Non-Overlapping Peak Set
#-------------------------------------------------------------------------------------------------
message("INFO : Make Non-Overlapping Peak Set")
df <- data.frame(
  samples = gsub("\\_summits.bed","",list.files(dirPeaks, pattern = "\\_summits.bed", full.names = FALSE)),
  groups = "scATAC",
  summits = list.files(dirPeaks, pattern = "\\_summits.bed", full.names = TRUE)
  )

unionPeaks <- extendedPeakSet(
    df = df,
    BSgenome = genome, 
    extend = 250,
    blacklist = "/home/ye/Work/BioAligment/SNP/Shi/blacklist/hg38_blacklist.bed",
    nSummits = 200000
  )
unionPeaks <- unionPeaks[seqnames(unionPeaks) %in% paste0("chr",c(1:22,"X"))]
unionPeaks <- keepSeqlevels(unionPeaks, paste0("chr",c(1:22,"X")))

saveRDS(unionPeaks,file.path(outDir,"unionPeaks.rds"))


message("INFO : Create Counts list")
countsPeaksList <- lapply(seq_along(bwFiles), function(i){
  message(sprintf("%s of %s", i, length(bwFiles)))
  gc()
  countInsertions(unionPeaks, import(bwFiles[i]),name=Names[i],by = "score")
 })

message("INFO : CountsMatrix")
mat <- lapply(countsPeaksList, function(x) x[[1]]) %>% Reduce("cbind",.)
frip <- lapply(countsPeaksList, function(x) x[[2]]) %>% unlist
total <- lapply(countsPeaksList, function(x) x[[3]]) %>% unlist

message("INFO : Create SummarizedExperiment Object")
se <- SummarizedExperiment(
  assays = SimpleList(counts = mat),
  rowRanges = unionPeaks
  )

Sample=unlist(lapply(colnames(se),function(name)return(stringr::str_split(name,"-")[[1]][1])))
rownames(se) <- paste(seqnames(se),start(se),end(se),sep="_")
colData(se)$FRIP <- frip
colData(se)$uniqueFrags <- total / 2
colData(se)$Sample=Sample

message("INFO : Save ATAC-Summarized-Experiment")
saveRDS(se,file.path(outDir,"scATAC-Summarized-Experiment.rds"))
