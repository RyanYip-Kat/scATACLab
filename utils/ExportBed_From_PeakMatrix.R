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
    useMatrix = "PeakMatrix",
    groupBy = args$groupby,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
peaks=getPeakSet(projHeme)
if(!is.null(peaks)){
        rowRanges(se)=peaks
        # mcols(peaks)
        rowData(se)=cbind(rowData(se),mcols(peaks))
        rgs=as.data.frame(ranges(peaks))
        sqs=as.data.frame(seqnames(peaks))
        peakset=cbind(sqs,rgs)
        peak_set=paste(peakset$value,peakset$start,peakset$end,sep="_")
        rownames(se)=peak_set
}

all_Peaks=rowData(se)
out <- data.frame(
              chr = all_Peaks$seqnames,
              start = c(as.integer(all_Peaks$start) - 1, as.integer(all_Peaks$start) - 1),
              end = c(as.integer(all_Peaks$end) - 1, as.integer(all_Peaks$end) - 1)
              )
path=file.path(args$outdir,"all_Peaks.bed")
cat(sprintf("### Write all Peaks into : %s \n",path))
readr::write_tsv(x = out,
          append = TRUE,
          path =path ,
          col_names = FALSE)


saveRDS(se,file.path(args$outdir,paste0(args$groupby,"_Differential_Summary.rds")))
markerList <- getMarkers(se, cutOff = "FDR <= 0.01 & Log2FC >= 1.25",returnGR=TRUE)
for(name in names(markerList)){
	path=file.path(args$outdir, paste0(str_replace(name," ","-"), ".bed"))
	fragmentsj <- markerList[[name]]
	cat(sprintf("### Write : %s Different peaks into : %s \n",name,path))
	#message(sprintf("### Write : %s Different peaks into : %s bed file\n",name,path))
	if(length(fragmentsj) > 0){
	    out <- data.frame(
	      chr = c(seqnames(fragmentsj), seqnames(fragmentsj)),
	      start = c(as.integer(start(fragmentsj) - 1), as.integer(end(fragmentsj) - 1)),
	      end = c(as.integer(start(fragmentsj)), as.integer(end(fragmentsj)))
	      ) %>% readr::write_tsv(
          x = .,
          append = TRUE,
          path = path,
          col_names = FALSE)
	}
}
