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

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 


print("### Loading ArrowFiles")
projHeme5<-loadArchRProject(args$project)

message("INFO :  Pseudo-time heatmaps ...")
trajMM  <- getTrajectory(ArchRProj = projHeme5, name = args$trajectory, useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))

trajGSM <- getTrajectory(ArchRProj = projHeme5, name = args$trajectory, useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))

trajPM  <- getTrajectory(ArchRProj = projHeme5, name = args$trajectory, useMatrix = "PeakMatrix", log2Norm = TRUE)
p3 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))

MyplotPDF(p1,p2,p3, name = paste0("Traj-Heatmap-",args$trajectory,".pdf"), outpath=args$outdir, addDOC = FALSE, width = 12, height = 16)

message("INFO : correlateTrajectories [ Motif and GeneScoreMatrix...")
corGSM_MM <- correlateTrajectories(trajGSM, trajMM,varCutOff1=0.6,varCutOff2=0.6,corCutOff=0.45)
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]

trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)


ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
MyplotPDF(ht1+ht2, name = paste0("TrajectoryHeatmapCorrelated-",args$trajectory,".pdf"), outpath=args$outdir, addDOC = FALSE, width = 16, height = 24)
