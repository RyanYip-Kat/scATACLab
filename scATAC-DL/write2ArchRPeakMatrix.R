library(argparse)
library(ArchR)
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
#########################  whether subset 
if(!is.null(args$column) & !is.null(args$subset)){
	message("INFO : Subset")
	metadata=as.data.frame(getCellColData(projHeme))
	cells=getCellNames(projHeme)
	stopifnot(args$column%in%colnames(metadata))
	target=metadata[[args$column]]
	
	stopifnot(args$subset%in%unique(target))
	idxPass= which(target%in%args$subset)
	cellPass=cells[idxPass]
	projHeme=subsetCells(projHeme,cellNames=cellPass)
}


message("INFO : get PeakMatrix")
se=getMatrixFromProject(projHeme,useMatrix="PeakMatrix")
peaks=getPeakSet(projHeme)
rowRanges(se)=peaks
rowData(se)=cbind(rowData(se),mcols(peaks))
peakNames=paste(seqnames(peaks),start(peaks),end(peaks),sep="_")
rownames(se)=peakNames

if(!is.null(args$peaktype)){
        se=subset(se,peakType%in%args$peaktype)
}


mat=assay(se)
message("INFO : Write PeakMatrix into 10X matrix")
filename=file.path(args$outdir,"ArchRPeakMatrix")
write10xCounts(x =mat, path=filename)


