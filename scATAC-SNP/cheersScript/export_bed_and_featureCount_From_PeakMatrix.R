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

parser$add_argument("--groupby",
                    type="character",
		    default="Clusters",
                    help="the column in cellcoldata in archr as group")

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
all_Peaks=rowData(se)
#if(!is.null(peaks)){
#        rowRanges(se)=peaks
#        # mcols(peaks)
#        rowData(se)=cbind(rowData(se),mcols(peaks))
#        rgs=as.data.frame(ranges(peaks))
#        sqs=as.data.frame(seqnames(peaks))
#        peakset=cbind(sqs,rgs)
#        peak_set=paste(peakset$value,peakset$start,peakset$end,sep="_")
#        rownames(se)=peak_set
#}
#out <- data.frame(
#              chr = c(all_Peaks$seqnames,all_Peaks$seqnames),
#              start = c(as.integer(all_Peaks$start) - 1, as.integer(all_Peaks$end) - 1),
#              end = c(as.integer(all_Peaks$start), as.integer(all_Peaks$end))
#              )

print("INFO : getPeakSets")
metadata=as.data.frame(getCellColData(projHeme))
cells=getCellNames(projHeme)
peaks=getPeakSet(projHeme)
Barcodes=rownames(metadata)

out=cbind(data.frame(chr=seqnames(peaks)),as.data.frame(ranges(peaks)))
out$dissToTSS=mcols(peaks)[["distToTSS"]]
out$Reproducibility=mcols(peaks)[["Reproducibility"]]
out$idx=paste0("Pbmc.",mcols(peaks)[["idx"]])
out$strand=as.character(strand(peaks))
out=out[,c("chr","start","end","idx","dissToTSS","strand","Reproducibility")]

path=file.path(args$outdir,"AllPeaks.bed")
cat(sprintf("INFO: Write all Peaks into : %s \n",path))
readr::write_tsv(x = out,
          append = FALSE,
          path =path ,
          col_names = FALSE)


saveRDS(se,file.path(args$outdir,paste0(args$groupby,"_Differential_Summary.rds")))
markerList <- getMarkers(se, cutOff = "FDR < 0.1 & Log2FC >0.25",returnGR=TRUE)
for(name in names(markerList)){
	path=file.path(args$outdir, paste0(str_replace(name," ","-"), ".bed"))
	fragmentsj <- markerList[[name]]
	cat(sprintf("### Write : %s Different peaks into : %s \n",name,path))
	#message(sprintf("### Write : %s Different peaks into : %s bed file\n",name,path))
	if(length(fragmentsj) > 0){
	    #out <- data.frame(
	    #  chr = c(seqnames(fragmentsj), seqnames(fragmentsj)),
	    #  start = c(as.integer(start(fragmentsj) - 1), as.integer(end(fragmentsj) - 1)),
	    #  end = c(as.integer(start(fragmentsj)), as.integer(end(fragmentsj)))
	    #  ) 
	  out=cbind(data.frame(chr=seqnames(fragmentsj)),
		    as.data.frame(ranges(fragmentsj))
		    )
	  out=out[,c("chr","start","end")]
	  readr::write_tsv(
		    x = out,
                    append = FALSE,
                    path = path,
                    col_names = FALSE)
	}

       #groups=unique(metadata$label_main)
       #cells=row.names(subset(metadata,label_main%in%groups[1]))
       cells=Barcodes[metadata[[args$groupby]]==name]
       cells=unlist(lapply(cells,function(x)str_split(x,"#")[[1]][2]))
       cells=unlist(lapply(cells,function(x)str_split(x,"-")[[1]][1]))

       filename=file.path(args$outdir, paste0(str_replace(name," ","-"), "_barcode.fa"))
       for(cell in cells){
	       cat(paste0(">",cell),sep="\n",file=filename,append=TRUE)
	       cat(cell,sep="\n",file=filename,append=TRUE)
       }
}
######################
message("### Get featureCount by group")
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

message("### Get PeakMatrix")
se=getMatrixFromProject(projHeme,useMatrix="PeakMatrix")
rowData(se)=peaks
all.out=cbind(data.frame(chr=seqnames(peaks)),as.data.frame(ranges(peaks)))
peak=paste(all.out$chr,all.out$start,all.out$end,sep="-")
rownames(se)=peak
mat=assay(se)
if(args$binarize){
        message("Making PseudoBulk...")
        mat@x[mat@x > 0] <- 1 #binarize
}
cluster=colData(se)[[args$groupby]]
clusterSums <- groupSums(mat = mat, groups = cluster, sparse = TRUE) #Group Sums

for(name in names(markerList)){
        fragmentsj <- markerList[[name]]
	out=cbind(data.frame(chr=seqnames(fragmentsj)),
                    as.data.frame(ranges(fragmentsj)))
	peak=paste(out$chr,out$start,out$end,sep="-")
	clustesum=clusterSums[peak,name]
	out$count=clustesum
	out=out[,c("chr","start","end","count")]
	path=file.path(args$outdir, paste0(str_replace(name," ","-"), "_ReadsInPeaks.txt"))
	write.table(out,path,quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)

        clustesum=clusterSums[,name]
	out=cbind(all.out,data.frame(count=clustesum))
	out=out[,c("chr","start","end","count")]
	path=file.path(args$outdir, paste0("total_",str_replace(name," ","-"), "_ReadsInPeaks.txt"))
	write.table(out,path,quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)	
}


