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
		    help="the project path  of ArchR")



parser$add_argument("--useMatrix",
                    type="character",
                    default="GeneScoreMatrix",
		    help="the matrix use in FindMarkers:PeakMatrix,MotifMatrix")

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--groupby",
                    type="character",
		    default="Clusters",
                    help="the column in cellcoldata in archr as group")


parser$add_argument("--column",
                    type="character",
                    default=NULL,
		    help="which column to subset")

parser$add_argument("--subset",
                    type="character",
		    nargs="+",
                    default=NULL)

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
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



print("# Get Marker Features in GeneMatrix slot")
se <- getMarkerFeatures(
    ArchRProj = projHeme,
    useMatrix = args$useMatrix,
    groupBy = args$groupby,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
if(args$useMatrix=="GeneScoreMatrix"){
	features=getFeatures(projHeme,useMatrix ="GeneScoreMatrix")
        rownames(se)=features
}else if(args$useMatrix=="PeakMatrix"){
	peaks=getPeakSet(projHeme)
        rowRanges(se)=peaks
        # mcols(peaks)
        rowData(se)=cbind(rowData(se),mcols(peaks))
        rgs=as.data.frame(ranges(peaks))
        sqs=as.data.frame(seqnames(peaks))
        peakset=cbind(sqs,rgs)
        peak_set=paste(peakset$value,peakset$start,peakset$end,sep="_")
        rownames(se)=peak_set
}else{
	motif=getMatrixFromProject(projHeme,useMatrix="MotifMatrix")
	rownames(se)=rownames(motif)
}


saveRDS(se,file.path(args$outdir,paste0(args$useMatrix,"_Differential_Summary.rds")))
markerList <- getMarkers(se, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.table(markerList,file.path(args$outdir,"MarkerList.csv"),sep=",",quote=F,row.names=F)

