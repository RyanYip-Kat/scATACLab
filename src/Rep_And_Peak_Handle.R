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


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)
projHeme2 <- addGroupCoverages(ArchRProj = projHeme, groupBy =args$groupby,force = TRUE)

pathToMacs2 <- findMacs2()
print(paste0("The Macs2 path is :",pathToMacs2))

print("# Add ReproduciblePeakSet")
projHeme2 <- addReproduciblePeakSet(
    ArchRProj = projHeme2,
    additionalParams = "--bdg --nomodel --nolambda",
    extendSummits=250,     # control peaks width
    groupBy = args$groupby,
    pathToMacs2 = pathToMacs2,
    force=TRUE
)

print("# Add Peaks Matrix")
projHeme2 <- addPeakMatrix(projHeme2)
print(getAvailableMatrices(projHeme2))
if(is.null(getImputeWeights(projHeme2))){
        projHeme2=addImputeWeights(projHeme2)
}

saveArchRProject(ArchRProj = projHeme2, outputDirectory = file.path(args$outdir,"Save-ProjHeme-PeakSetMacs2"), load = FALSE)

