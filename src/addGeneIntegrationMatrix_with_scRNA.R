library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(viridis)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
		    default=NULL,
                    help="the path of project saved")


parser$add_argument("--scRNA",
                    type="character",
                    default=NULL,
		    help="the scRNA object : SE or Seurat")


parser$add_argument("--groupATAC",
                    type="character",
                    default="Clusters",
                    help="the column in ATAC colData")


parser$add_argument("--groupRNA",
                    type="character",
                    default="seurat_clusters",
                    help="the column in RNA colData or @meta.data")


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

print("### Loading ArrowFiles")
projHeme5<-loadArchRProject(args$project)
seRNA<-readRDS(args$scRNA)

print("### addGeneIntegrationMatrix")
projHeme5 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme5,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    groupRNA = args$groupRNA,
    groupATAC=args$groupATAC,
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force=TRUE,
)

saveArchRProject(ArchRProj = projHeme5, outputDirectory = file.path(args$outdir,"Save-ProjHeme-Align-scRNA"), load =TRUE)

