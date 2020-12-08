library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
library(Cairo)
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the project path  of ArchR")


parser$add_argument("--markers",
		    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--groupby",
                    type="character",
                    default="Clusters")

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

markerGenes=args$markers
print(paste0("# The plot size is : [ ",width,height," ]"))

### add own region
#df=data.frame(chrom="chr11",start=c(100157509,100187509),end=100207509,100257510)
#gr3 <- as(df, "GRanges")
p <- plotBrowserTrack(ArchRProj = projHeme,
			      groupBy =args$groupby,
                              geneSymbol = markerGenes,
			      #region=gr3,
                              upstream = 50000,
                              downstream = 50000)
MyplotPDF(plotList = p,
		name = "Plot-Tracks-Marker-Genes.pdf",
                outpath = args$outdir,
                addDOC = FALSE, width = width, height = height)


