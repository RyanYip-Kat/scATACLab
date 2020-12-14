library(argparse)
library(stringr)
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

parser$add_argument("--normal",
                    type="character",
                    default="HC",
		    help="the column in Status as normal group")

parser$add_argument("--effect",
                    type="character",
                    default="VKH",
		    help="the column in Status as effect group")


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
cells=getCellNames(projHeme)

print("INFO : Subset from normal and effect")
target=metadata[["Status"]]
idxPass= which(target%in%c(args$normal,args$effect))
cellPass=cells[idxPass]
projHeme=subsetCells(projHeme,cellNames=cellPass)

#if(!is.null(args$column) & !is.null(args$subset)){
#	stopifnot(args$column%in%colnames(metadata))
#	target=metadata[[args$column]]
#	
#	stopifnot(args$subset%in%unique(target))
#	idxPass= which(target%in%args$subset)
#	cellPass=cells[idxPass]
#	projHeme=subsetCells(projHeme,cellNames=cellPass)
#}


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
mat=assay(se)
saveRDS(se,file.path(args$outdir,"scATAC-GeneSummarized-Experiment.rds"))

#rownames(mat)=paste(as.character(rowData(se)[["GroupReplicate"]]),1:nrow(mat),sep="#")
metadata=as.data.frame(getCellColData(projHeme))
cells=getCellNames(projHeme)
Barcodes=rownames(metadata)
stopifnot(args$column%in%colnames(metadata))
target=metadata[[args$column]]


message("INFO :Get Gene featureCount by group")
counts_list=list()
for(sub in args$subset){
        idxPass= which(target%in%sub)
        #cellPass=cells[idxPass]
	status=as.character(metadata[idxPass,][["Status"]])
	cluster=as.character(metadata[idxPass,][["Sample"]])
	group=paste(cluster,status,sep="_")
	mm=mat[,idxPass]
	clusterSums <- groupSums(mat = mm, groups = group, sparse = TRUE) #Group Sums
	colnames(clusterSums)=paste0(sub,"_",colnames(clusterSums))
	counts_list[[sub]]=clusterSums
}

#counts_list=list()
#for(sub in c("DC","BC","NK","TC")){
#        idxPass= which(target%in%sub)
#        cellPass=cells[idxPass]
#        status=as.character(metadata[idxPass,][["Status"]])
#        cluster=as.character(metadata[idxPass,][["Sample"]])
#        group=paste(cluster,status,sep="_")
#        mm=mat[,idxPass]
#        clusterSums <- groupSums(mat = mm, groups = group, sparse = TRUE) #Group Sums
#        colnames(clusterSums)=paste0(sub,"_",colnames(clusterSums))
#        counts_list[[sub]]=clusterSums
#}

message("INFO : Save GenesCounts")
counts=do.call(cbind,counts_list)
path=file.path(args$outdir,"scATAC-GeneScoreCounts.txt")
write.table(counts,path,sep="\t",quote=FALSE)

