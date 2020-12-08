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
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

parser$add_argument("--groupby",
                    type="character",
                    default="Clusters")

parser$add_argument("--name",
                    type="character",
                    default="CEpi-0h-3h")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("### Loading ArrowFiles")
projHeme5<-loadArchRProject(args$project)

print("### addTrajectory")
projHeme5 <- addTrajectory(
    ArchRProj = projHeme5, 
    name = args$name, 
    groupBy = args$groupby,
    trajectory = c("C1","C2","C3"), 
    embedding = "UMAP", 
    force = TRUE
)

saveArchRProject(ArchRProj = projHeme5, outputDirectory = file.path(args$outdir,"Save-ProjHeme-Trajectory"), load = TRUE)
