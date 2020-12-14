library(argparse)
library(stringr)
library(Matrix)
library(ggplot2)
library(edgeR)
library(matrixStats)
library(SummarizedExperiment)
library(readr)
library(GenomicRanges)
library(magrittr)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--se",
                    type="character",
                    default=NULL,
		    help="SummarizedExperiment object")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--binarize",
                    action="store_true")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
		    help="which column to subset")

parser$add_argument("--subset",
                    type="character",
		    nargs="+",
                    default=NULL)

args <- parser$parse_args()

outDir=args$outdir
if(!dir.exists(outDir)){
        dir.create(outDir,recursive=TRUE)
}

message("INFO : Loading dataset")
se=readRDS(args$se)
meta=as.data.frame(colData(se))
if(!is.null(args$column) & !is.null(args$subset)){
	col=meta[[args$column]]
	idx=col%in%args$subset
	se=se[,idx]
}

######################
message("INFO :Get featureCount by group")
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

message("INFO : Get PeakMatrix")
mat=assay(se)
if(args$binarize){
        message("Making PseudoBulk...")
        mat@x[mat@x > 0] <- 1 #binarize
}

Sample=colData(se)$Sample
counts=groupSums(mat,Sample,sparse=TRUE)

message("INFO : Save featureCounts")
path=file.path(outDir,"scATAC-featureCounts.txt")
write.table(counts,path,sep="\t",quote=FALSE)

out=data.frame("chr"=seqnames(se),as.data.frame(ranges(se)))
out=out[,c("chr","start","end")]
write.table(out,file.path(outDir,"AllPeaks.bed"),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

