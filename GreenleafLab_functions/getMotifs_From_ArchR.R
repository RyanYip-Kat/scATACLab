library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
library(edgeR)
library(matrixStats)
library(SummarizedExperiment)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

parser$add_argument("--nTop",
                    type="integer",
                    default=50000,
                    help="number of top var peaks")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print(paste0("### Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("### Loading Project")
projHeme<-loadArchRProject(args$project)

avaliable_matrixs=getAvailableMatrices(projHeme)
se=getMatrixFromProject(projHeme,useMatrix="MotifMatrix")


print("### Add ReducedDims")
ReducedDims_Names=c("IterativeLSI","Harmony")
for(name in ReducedDims_Names){
	slot=getReducedDims(projHeme,reducedDims=name,returnMatrix=F)
	metadata(se)[[name]]=slot
}

print(names(metadata(se)))
print("### Add Embedding")
mat=as.matrix(getEmbedding(ArchRProj =projHeme, embedding ="UMAP", returnDF = TRUE))
#Add UMAP coordinates to column data in summarized experiment
colData(se)$UMAP1 <-mat[,1]
colData(se)$UMAP2 <-mat[,2]

mat=as.matrix(getEmbedding(ArchRProj =projHeme, embedding ="TSNE", returnDF = TRUE))
colData(se)$TSNE1 <-mat[,1]
colData(se)$TSNE2 <-mat[,2]


###################
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
#Make Pseudo Bulk Library
message("Making PseudoBulk...")
mat=assay(se)
mat@x[mat@x > 0] <- 1 #binarize
cluster=colData(se)$Clusters
clusterSums <- groupSums(mat = mat, groups = cluster, sparse = TRUE) #Group Sums
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), args$nTop) #Top variable peaks
metadata(se)$variablePeaks <- varPeaks


print("### Save Matrix")
saveRDS(se,file.path(args$outdir,"MotifMatrix_Summarized-Experiment.rds"))
