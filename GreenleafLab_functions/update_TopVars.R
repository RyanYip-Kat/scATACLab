library(argparse)
library(stringr)
library(Matrix)
library(SummarizedExperiment)
library(edgeR)
library(matrixStats)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--se",
                    type="character",
                    default="SummarizedExperiment object")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--nTop",
                    type="integer",
                    default=50000 ,
                    help="the top peaks to be used")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

###############
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
###############
print("### Loading dataset")
se=readRDS(args$se)

#Make Pseudo Bulk Library
message("Making PseudoBulk...")
mat=assay(se)
mat@x[mat@x > 0] <- 1 #binarize
cluster=colData(se)$Clusters
clusterSums <- groupSums(mat = mat, groups = cluster, sparse = TRUE) #Group Sums
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), args$nTop) #Top variable peaks
metadata(se)$variablePeaks <- varPeaks

print("### Save...")
saveRDS(se,args$se)
