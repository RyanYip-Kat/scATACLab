library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
library(Seurat)
library(stringr)
library(argparse)

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


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

message("INFO : loading dataset ...")
se=readRDS(args$se)
######################
message("### Get featureCount by group")
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

message("INFO : get counts")
mat=assay(se)
if(args$binarize){
        message("Making PseudoBulk...")
        mat@x[mat@x > 0] <- 1 #binarize
}
cluster=colData(se)[["Sample"]]
clusterSums <- groupSums(mat = mat, groups = cluster, sparse = TRUE) #Group Sums

all.out=data.frame("chr"=seqnames(se),as.data.frame(ranges(se)))
all.Count=rowSums(clusterSums)
for(name in unique(cluster)){
	cat(sprintf("Sample : [ %s ] \n",name))
        clustesum=clusterSums[,name]
	out=cbind(all.out,data.frame(count=clustesum))
	out=out[,c("chr","start","end","count")]
	path=file.path(args$outdir, paste0(str_replace(name," ","-"), "_ReadsInPeaks.txt"))
	write.table(out,path,quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)	
}

out=cbind(all.out,data.frame(count=all.Count))
out=out[,c("chr","start","end","count")]
path=file.path(args$outdir, "All_ReadsInPeaks.txt")
write.table(out,path,quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)


