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

parser$add_argument("--trajectory",
		    nargs="+",
		    type="character",
		    default=NULL,
		    help="subset name in groupby")

parser$add_argument("--name",
                    type="character",
                    default="test",
		    help="Trajectory Name")


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
    trajectory = args$trajectory, 
    embedding = "UMAP", 
    reducedDims = "IterativeLSI",
    force = TRUE
)

outName=paste0("Save-ProjHeme-",args$name,"-Trajectory")
saveArchRProject(ArchRProj = projHeme5, outputDirectory = file.path(args$outdir,outName), load = TRUE)
