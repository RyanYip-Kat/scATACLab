library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
		    default=NULL,
                    help="the path of project saved")

parser$add_argument("--group",
                    type="character",
                    default="Clusters",
		    help="the column in cellcol to be exported")

parser$add_argument("--trajectory",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset of column to be extracted")

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--name",
                    type="character",
                    default="name for trajectoty slot")

parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading data")
projHeme=loadArchRProject(args$project)

print("# add Trajectory")
trajectory=args$trajectory
projHeme <- addTrajectory(
    ArchRProj = projHeme,
    name = args$name,
    groupBy = args$group,
    trajectory = trajectory,
    embedding = "UMAP",
    force = TRUE
)

print("# Save ...")
saveArchRProject(ArchRProj = projHeme, outputDirectory =file.path(args$outdir,"Save-ProjHeme"), load = TRUE)

p <- plotTrajectory(projHeme, trajectory =args$name, colorBy = "cellColData", name = args$name)
plotPDF(p, name = paste0("Plot-",args$name,"-Traj-UMAP.pdf"), ArchRProj = projHeme, addDOC = FALSE, width =12, height = 12)

print("# addImputeWeights")
projHeme <- addImputeWeights(projHeme)

print("# Pstime heatmap")
trajMM  <- getTrajectory(ArchRProj = projHeme, name = args$name, useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
trajGSM <- getTrajectory(ArchRProj = projHeme, name = args$name, useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
trajPM  <- getTrajectory(ArchRProj = projHeme5, name = args$name, useMatrix = "PeakMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))

plotPDF(p1, p2, p3, name = paste0("Plot-",args$name,"-Traj-Heatmaps.pdf"), ArchRProj = projHeme, addDOC = FALSE, width = 12, height = 16)

