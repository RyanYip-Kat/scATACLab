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
                    default="the path of project saved,the project intergrated with scRNA,had GeneIntegrationMatrix in AvailableMatrices")


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



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("### Loading ArrowFiles")
projHeme5<-loadArchRProject(args$project)

print("### addCoAccessibility")
projHeme5 <- addCoAccessibility(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
    ArchRProj = projHeme5,
    corCutOff = 0.3,
    resolution = 1,
    returnLoops = TRUE
)


saveRDS(cA,file.path(args$outdir,"CoAccessibility.rds"))
saveArchRProject(ArchRProj = projHeme5, outputDirectory = file.path(args$outdir,"Save-ProjHeme-CoAccessibility"), load = TRUE)
