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


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
projHeme=loadArchRProject(args$project)


width=args$width
height=args$height

print(paste0("# The plot size is : [ ",width,height," ]"))
p1 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = args$groupby, embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = "Sample", embedding = "UMAP")


MyplotPDF(p1,p2, name = paste0("Plot-UMAP-Sample-",args$groupby,".pdf"), outpath=args$outdir, addDOC = FALSE, width = width, height = height)

p1 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name =args$groupby, embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = projHeme, colorBy = "cellColData", name = "Sample", embedding = "TSNE")

MyplotPDF(p1,p2, name = paste0("Plot-TSNE-Sample-",args$groupby,".pdf"), outpath=args$outdir, addDOC = FALSE, width = width, height = height)

