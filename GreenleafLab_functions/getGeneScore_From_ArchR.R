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


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print(paste0("### Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("### Loading Project")
projHeme<-loadArchRProject(args$project)

avaliable_matrixs=getAvailableMatrices(projHeme)
stopifnot("GeneScoreMatrix"%in%avaliable_matrixs)
se=getMatrixFromProject(projHeme,useMatrix="GeneScoreMatrix")

features=getFeatures(projHeme,useMatrix ="GeneScoreMatrix")
rownames(se)=features
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


print("### Save  Matrix")
saveRDS(se,file.path(args$outdir,"GeneScoreMatrix_Summarized-Experiment.rds"))
