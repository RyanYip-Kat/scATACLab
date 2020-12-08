library(argparse)
library(stringr)
library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(viridis)

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
                    default=4,
                    help="number of threads for command")


parser$add_argument("--groupRNA",
                    type="character",
                    default="label.fine")

parser$add_argument("--scrna",
                    type="character",
                    default=NULL)

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)

print("# Add GeneIntegrationMatrix")
print(paste0("# Finetune with scRNA: ",args$scrna))
seRNA=readRDS(args$scrna)
projHeme2 <- addGeneIntegrationMatrix(
    ArchRProj = projHeme,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = args$groupRNA,
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(projHeme2$Clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
print(cbind(preClust, rownames(cM))) #Assignments

print(unique(unique(projHeme2$predictedGroup_Un)))


saveArchRProject(ArchRProj = projHeme2, outputDirectory =file.path(args$outdir,"Save-ProjHeme-Prediction_Un"), load =FALSE)

