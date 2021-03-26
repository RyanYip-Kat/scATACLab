library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)

source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Trajectory plot')
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

parser$add_argument("--trajectory",
		    type="character",
		    default=NULL,
		    help="trajectory name")

parser$add_argument("--markers",
		    nargs="+",
		    type="character",
		    default=NULL,
		    help="genes to plot in trajectory")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 


print("### Loading ArrowFiles")
projHeme5<-loadArchRProject(args$project)

p <- plotTrajectory(projHeme5, trajectory = args$trajectory, colorBy = "cellColData", name = args$trajectory)
MyplotPDF(p, name = paste0("Trajectory-",args$trajectory,".pdf"), outpath=args$outdir, addDOC = FALSE, width =12, height = 12)

if(!is.null(args$markers)){
	if(is.null(getImputeWeights(projHeme5))){
		projHeme=addImputeWeights(projHeme5)
	}
	for(gene in args$markers){
		p <- plotTrajectory(projHeme5, trajectory = args$trajectory, colorBy = "GeneScoreMatrix", name = gene, continuousSet = "horizonExtra")
		MyplotPDF(p, name = paste0("TrajectoryGene-",gene,".pdf"), outpath=args$outdir, addDOC = FALSE, width = 16, height = 12)
	}
}

message("INFO :  Pseudo-time heatmaps ...")
trajMM  <- getTrajectory(ArchRProj = projHeme5, name = args$trajectory, useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))

trajGSM <- getTrajectory(ArchRProj = projHeme5, name = args$trajectory, useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))

trajPM  <- getTrajectory(ArchRProj = projHeme5, name = args$trajectory, useMatrix = "PeakMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))

MyplotPDF(p1,p2,p3, name = paste0("Traj-Heatmap-",args$trajectory,".pdf"), outpath=args$outdir, addDOC = FALSE, width = 12, height = 16)

