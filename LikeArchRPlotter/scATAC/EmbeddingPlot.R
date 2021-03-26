library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the project path  of ArchR")


parser$add_argument("--groupby",
                    type="character",
                    default="Clusters",
                    help="the project path  of ArchR")

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--width",
                    type="integer",
                    default=12,
                    help="the width of plot")

parser$add_argument("--height",
                    type="integer",
                    default=12,
                    help="the height of plot")

parser$add_argument("--column",
                    type="character",
                    default="Clusters",
                    help="the column in cellcol to be exported")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset of column to be extracted")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
projHeme=loadArchRProject(args$project)
metadata=as.data.frame(getCellColData(projHeme))
stopifnot(args$groupby%in%colnames(metadata))
cells=getCellNames(projHeme)
if(!is.null(args$column) & !is.null(args$subset)){
        stopifnot(args$column%in%colnames(metadata))
        target=metadata[[args$column]]

        stopifnot(args$subset%in%unique(target))
        idxPass= which(target%in%args$subset)
        cellPass=cells[idxPass]
        projHeme=subsetCells(projHeme,cellNames=cellPass)
}


width=args$width
height=args$height

p1 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = args$groupby, embedding = "UMAP")
#p2 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = "Sample", embedding = "UMAP")


MyplotPDF(p1, name = paste0("Plot-UMAP-Sample-",args$groupby,".pdf"), outpath=args$outdir, addDOC = FALSE, width = width, height = height)

p1 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name =args$groupby, embedding = "TSNE")
#p2 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = "Sample", embedding = "TSNE")

MyplotPDF(p1, name = paste0("Plot-TSNE-Sample-",args$groupby,".pdf"), outpath=args$outdir, addDOC = FALSE, width = width, height = height)

