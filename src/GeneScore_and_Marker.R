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

parser$add_argument("--markers",
		    nargs="+",
                    type="character",
                    default=NULL)

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)

print("# Get Marker Features")
markersGS <- getMarkerFeatures(
    ArchRProj = projHeme, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = args$groupby,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
out_path=getOutputDirectory(projHeme)
saveRDS(markersGS,file.path(out_path,paste0(args$groupby,"_","markersGS.rds")))


markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.table(markerList,file.path(args$outdir,"markerList.csv"),sep=",",quote=F,row.names=F)

markerGenes=args$markers
if(!is.null(markerGenes)){
	print(paste0("### Plot heatmap : ",length(markerGenes)," genes"))
}

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 16, height = 12, ArchRProj = projHeme, addDOC = FALSE)
