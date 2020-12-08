library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
library(ggplot2)
library(edgeR)
library(matrixStats)
library(DropletUtils)
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

print("### Get Score GeneMatrix")
gs=getMatrixFromProject(projHeme,useMatrix = "GeneScoreMatrix")
features=getFeatures(projHeme,useMatrix = "GeneScoreMatrix")
rownames(gs)=features
mat=assay(gs)


message("### Rename Cells")
cells=getCellNames(projHeme)
barcodes=unlist(lapply(cells,function(cell)return(str_split(cell,"#")[[1]][2])))
barcodes=unlist(lapply(barcodes,function(cell)return(str_split(cell,"-")[[1]][1])))
idents=as.integer(as.factor(unlist(lapply(cells,function(cell)return(str_split(cell,"#")[[1]][1])))))
barcodes=paste(barcodes,idents,sep="-")
colnames(mat)=barcodes

message("### Update ArchR Project")
projHeme$Barcode=barcodes
path=getOutputDirectory(projHeme)
saveArchRProject(ArchRProj = projHeme, outputDirectory =path, load = TRUE)

message("### Write Gene Score Matrix into 10X matrix")
filename=file.path(args$outdir,"genescore_feature_bc_matrix.h5")
print(paste0("### Write counts matrix into : ",filename))
write10xCounts(x =mat, path=filename,type="HDF5")
filename=file.path(args$outdir,"GeneScoreMatrix")
print(paste0("### Write counts matrix into : ",filename))
write10xCounts(x =mat, path=filename)


