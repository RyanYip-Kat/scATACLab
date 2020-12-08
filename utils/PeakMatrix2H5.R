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


parser$add_argument("--binarize",
                    action="store_true")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
		    help="which column to subset")

parser$add_argument("--subset",
                    type="character",
		    nargs="+",
                    default=NULL)

parser$add_argument("--peaktype",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the peakType in PeakMatrix mcols for subset(Distal,Exonic,Intronic,Promoter)")

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



se=getMatrixFromProject(projHeme,useMatrix="PeakMatrix")
peaks=getPeakSet(projHeme)
if(!is.null(peaks)){
        rowRanges(se)=peaks
        # mcols(peaks)
        rowData(se)=cbind(rowData(se),mcols(peaks))
        rgs=as.data.frame(ranges(peaks))
        sqs=as.data.frame(seqnames(peaks))
        peakset=cbind(sqs,rgs)
	mcols(peaks)=cbind(mcols(peaks),peakset)
	peak_set=paste(peaks$value,peaks$start,peaks$end,sep="-")
	peak_set=paste0(peak_set,"(",peaks$nearestGene,")")
        rownames(se)=peak_set
}

######################
#message("### Get featureCount by group")
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

if(!is.null(args$peaktype)){
        se=subset(se,peakType%in%args$peaktype)
}


mat=assay(se)
if(args$binarize){
        message("### Making PseudoBulk...")
        mat@x[mat@x > 0] <- 1 #binarize
}
#cluster=colData(se)[[args$groupby]]
#clusterSums <- groupSums(mat = mat, groups = cluster, sparse = TRUE) #Group Sums
message("### Rename PeakMatrix")
cells=getCellNames(projHeme)
barcodes=unlist(lapply(cells,function(cell)return(str_split(cell,"#")[[1]][2])))
barcodes=unlist(lapply(barcodes,function(cell)return(str_split(cell,"-")[[1]][1])))
idents=as.integer(as.factor(unlist(lapply(cells,function(cell)return(str_split(cell,"#")[[1]][1])))))
barcodes=paste(barcodes,idents,sep="-")
colnames(mat)=barcodes

message("### Write PeakMatrix into 10X matrix")
filename=file.path(args$outdir,"filtered_peak_bc_matrix.h5")
print(paste0("### Write counts matrix into : ",filename))
write10xCounts(x =mat, path=filename,type="HDF5")
filename=file.path(args$outdir,"matrix")
print(paste0("### Write counts matrix into : ",filename))
write10xCounts(x =mat, path=filename)


