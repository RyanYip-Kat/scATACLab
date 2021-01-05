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
message("INFO : Done!")
