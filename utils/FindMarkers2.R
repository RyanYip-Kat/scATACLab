library(argparse)
library(stringr)
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


parser$add_argument("--useMatrix",
                    type="character",
                    default="GeneScoreMatrix",
		    help="the matrix use in FindMarkers")

parser$add_argument("--logfc",
		    type="double",
		    default=1.0,
		    help="logfc cutoff for getMarkers")

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--groupby",
                    type="character",
		    default="Clusters",
                    help="the column in cellcoldata in archr as group")

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
logfc=abs(args$logfc)
projHeme=loadArchRProject(args$project)

metadata=as.data.frame(getCellColData(projHeme))
stopifnot(args$groupby%in%colnames(metadata))
stopifnot(args$useMatrix%in%getAvailableMatrices(projHeme))

cells=getCellNames(projHeme)
if(!is.null(args$column) & !is.null(args$subset)){
        stopifnot(args$column%in%colnames(metadata))
        target=metadata[[args$column]]

        stopifnot(args$subset%in%unique(target))
        idxPass= which(target%in%args$subset)
        cellPass=cells[idxPass]
        projHeme=subsetCells(projHeme,cellNames=cellPass)
}

peakset=getPeakSet(projHeme)
message("INFO : Get Marker Features in GeneMatrix slot")
MF <- getMarkerFeatures(
    ArchRProj = projHeme,
    useMatrix = args$useMatrix,
    groupBy = args$groupby,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

message("INFO : Save ...")
saveRDS(MF,file.path(args$outdir,paste0("MF-",args$groupby,"-",args$useMatrix,".rds")))
markerUp <- getMarkers(MF, cutOff = paste0("Pval <= 0.05 & Log2FC >= ",logfc),returnGR=TRUE)
markerUpList<-list()
for(name in names(markerUp)){
	cat(sprintf("INFO : Up %s\n",name))
	gr=markerUp[[name]]
	grOverlap=subsetByOverlaps(peakset,gr)
	mcols(grOverlap)$Clusters=names(grOverlap)
	names(grOverlap)=paste(seqnames(grOverlap),start(grOverlap),end(grOverlap),sep="-")
	names(gr)=paste(seqnames(gr),start(gr),end(gr),sep="-")
	grDF=as.data.frame(gr)
	grOverlapDF=as.data.frame(grOverlap)
	grDF=grDF[rownames(grOverlapDF),]
	grOverlapDF=cbind(grOverlapDF,grDF[,c("Log2FC","FDR","MeanDiff")])
	grOverlapDF$Name=name
	grOverlapDF$peak=rownames(grOverlapDF)
	markerUpList[[name]]=grOverlapDF
}
m=do.call(rbind,markerUpList)
write.table(m,file.path(args$outdir,"MarkerUp.csv"),sep=",",quote=F,row.names=F)

markerDo <- getMarkers(MF, cutOff = paste0("Pval <= 0.05 & Log2FC <= ",-logfc),returnGR=T)
markerDoList<-list()
for(name in names(markerDo)){
        cat(sprintf("INFO : Up %s\n",name))
        gr=markerDo[[name]]
        grOverlap=subsetByOverlaps(peakset,gr)
        mcols(grOverlap)$Clusters=names(grOverlap)
        names(grOverlap)=paste(seqnames(grOverlap),start(grOverlap),end(grOverlap),sep="-")
        names(gr)=paste(seqnames(gr),start(gr),end(gr),sep="-")
        grDF=as.data.frame(gr)
        grOverlapDF=as.data.frame(grOverlap)
        grDF=grDF[rownames(grOverlapDF),]
        grOverlapDF=cbind(grOverlapDF,grDF[,c("Log2FC","FDR","MeanDiff")])
        grOverlapDF$Name=name
	grOverlapDF$peak=rownames(grOverlapDF)
        markerDoList[[name]]=grOverlapDF
}
m=do.call(rbind,markerDoList)
write.table(m,file.path(args$outdir,"MarkerDown.csv"),sep=",",quote=F,row.names=F)

heatmapMF <- plotMarkerHeatmap(
  seMarker = MF,
  cutOff = paste0("Pval <= 0.05 & Log2FC >= ",logfc),
  labelMarkers = NULL,
  labelRows=FALSE,
  transpose = TRUE
)
MyplotPDF(heatmapMF, name = paste0(args$useMatrix,"-Marker-Heatmap-Up"), width = 16, height =12, outpath =args$outdir, addDOC = FALSE)

heatmapMF <- plotMarkerHeatmap(
  seMarker = MF,
  cutOff = paste0("Pval <= 0.05 & Log2FC <= ",-logfc),
  labelMarkers = NULL,
  labelRows=FALSE,
  transpose = TRUE
)
MyplotPDF(heatmapMF, name = paste0(args$useMatrix,"-Marker-Heatmap-Do"), width = 16, height =12, outpath =args$outdir, addDOC = FALSE)

message("INFO : Done!")

