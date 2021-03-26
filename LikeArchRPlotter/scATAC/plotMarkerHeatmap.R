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


parser$add_argument("--markers",
		    nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--useMatrix",
                    type="character",
                    default="GeneScoreMatrix",
		    help="the matrix use in FindMarkers")

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--groupby",
                    type="character",
		    default="Clusters",
                    help="the column in cellcoldata in archr as group")

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
metadata=as.data.frame(getCellColData(projHeme))
stopifnot(args$groupby%in%colnames(metadata))
stopifnot(args$useMatrix%in%getAvailableMatrices(projHeme))

print("# Get Marker Features in GeneMatrix slot")
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme,
    useMatrix = args$useMatrix,
    groupBy = args$groupby,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerGenes=args$markers
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.table(markerList,file.path(args$outdir,"MarkerList.csv"),sep=",",quote=F,row.names=F)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
MyplotPDF(heatmapGS, name = paste0(args$useMatrix,"-Marker-Heatmap1"), width = width, height =height, outpath =args$outdir, addDOC = FALSE)
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = NULL,
  labelRows=TRUE,
  transpose = TRUE
)
MyplotPDF(heatmapGS, name = paste0(args$useMatrix,"-Marker-Heatmap2"), width = width, height =height, outpath =args$outdir, addDOC = FALSE)
#print("# Marker Peaks in Browser Tracks")
#p <- plotBrowserTrack(
#    ArchRProj = projHeme, 
#    groupBy = "Clusters2", 
#    geneSymbol = c("GATA1"),
#    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Erythroid"],
#    upstream = 50000,
#    downstream = 50000
#)
