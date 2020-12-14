library(argparse)
library(stringr)
library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
library(edgeR)
library(matrixStats)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the project path  of ArchR")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--column",
                    type="character",
                    default=NULL,
		    help="celltype column's name")

parser$add_argument("--subset",
                    type="character",
		    nargs="+",
                    default=NULL,
		    help="celltypes")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
projHeme=loadArchRProject(args$project)
metadata=as.data.frame(getCellColData(projHeme))
cells=getCellNames(projHeme)


if(!is.null(args$column) & !is.null(args$subset)){
	stopifnot(args$column%in%colnames(metadata))
	target=metadata[[args$column]]
	
	stopifnot(args$subset%in%unique(target))
	idxPass= which(target%in%args$subset)
	cellPass=cells[idxPass]
	projHeme=subsetCells(projHeme,cellNames=cellPass)
}


######################
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

message("INFO : Get Score GeneMatrix")
se=getMatrixFromProject(projHeme,useMatrix = "GeneScoreMatrix")
features=getFeatures(projHeme,useMatrix = "GeneScoreMatrix")
rownames(se)=features
colData(se)$celltype=colData(se)$label_fine
mat=assay(se)
saveRDS(se,file.path(args$outdir,"scATAC-GeneSummarized-Experiment.rds"))

#rownames(mat)=paste(as.character(rowData(se)[["GroupReplicate"]]),1:nrow(mat),sep="#")
metadata=as.data.frame(getCellColData(projHeme))
cells=getCellNames(projHeme)
Barcodes=rownames(metadata)

message("INFO : CreateSeurat")
seurat=CreateSeuratObject(mat,project="scATAC-GeneScore",meta.data=metadata)
message("INFO : Normalize")
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)

message("INFO : Save ")
saveRDS(seurat,file.path(args$outdir,"geneScore-seurat.rds"))
