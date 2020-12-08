library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(viridis)
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
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



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
projHeme5<-loadArchRProject(args$project)

print("# Motif")
if("Motif" %ni% names(projHeme5@peakAnnotation)){
    projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
    projHeme5 <- addBgdPeaks(projHeme5)
}

projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5,
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

print("# EncodeTFBS")
if("EncodeTFBS" %ni% names(projHeme5@peakAnnotation)){
    projHeme5 <- addArchRAnnotations(ArchRProj = projHeme5, collection = "EncodeTFBS")
}

plotVarDev <- getVarDeviations(projHeme5, plot = TRUE, name = "EncodeTFBSMatrix")
plotPDF(plotVarDev, name = "Variable-EncodeTFBS-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

print("# ATAC")
if("ATAC" %ni% names(projHeme5@peakAnnotation)){
    projHeme5 <- addArchRAnnotations(ArchRProj = projHeme5, collection = "ATAC")
}

projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5,
  peakAnnotation = "ATAC",
  force = TRUE
)

plotVarDev <- getVarDeviations(projHeme5, plot = TRUE, name = "ATACMatrix")
plotPDF(plotVarDev, name = "Variable-ATAC-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

saveArchRProject(ArchRProj = projHeme5, outputDirectory = file.path(args$outdir,"Save-ProjHeme-Motif"), load = TRUE)

