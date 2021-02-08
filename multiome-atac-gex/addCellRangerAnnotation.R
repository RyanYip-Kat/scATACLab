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
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--annotation",
                    type="character",
                    default=NULL,
		    help="*atac_peak_annotation.tsv")

parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



#print(paste0("Setting default genome  with :",args$genome))
#addArchRGenome(args$genome)

cat(sprintf("INFO : Setting threads  : %d\n",args$num_threads))
addArchRThreads(threads = args$num_threads) 

message("INFO :  Loading dataset ...")
proj<-loadArchRProject(args$project)
#if(is.null(getImputeWeights(proj))){
#        proj=addImputeWeights(proj)
#}

message("INFO : Loading Annotation ...")
keepChr=getSeqnames(proj)
AnnoDF=read.table(args$annotation,sep="\t",header=TRUE)
AnnoDF=subset(AnnoDF,str_detect(peak,paste(keepChr,collapse="|")))
message("INFO : make Annotation Granges ...")
peakList=list()
peakListGR=list()
for(i in 1:nrow(AnnoDF)){
	if(i%%100==0){
		cat(sprintf("INFO : [ %d of %d ]\n",i,nrow(AnnoDF)))
	}
	peak=AnnoDF$peak[i]
	pl=str_split(peak,"_|-")[[1]]
	x=data.frame("chr"=pl[1],"start"=pl[2],"end"=pl[3])
	peakList[[i]]=x
	#xx=cbind(x,AnnoDF[i,,drop=FALSE])
	#gr=makeGRangesFromDataFrame(xx,keep.extra.columns=TRUE)
	#peakListGR[[i]]=gr
}
peakDF=do.call(rbind,peakList)
AnnoDF=cbind(peakDF,AnnoDF)
AnnoGR=makeGRangesFromDataFrame(AnnoDF,keep.extra.columns=TRUE)
AnnoGRList=GRangesList("CRG"=AnnoGR)

message("INFO : addPeakSet ...")
proj=addPeakSet(proj,peakSet=AnnoGR,force=TRUE)

message("INFO : addPeakAnnotations ...")
proj=addPeakAnnotations(proj,regions=AnnoGRList,name="CellRanger")

message("INFO : Save ..")
saveArchRProject(ArchRProj = proj,overwrite=TRUE,outputDirectory =file.path(args$outdir,"Save-ProjHeme-CRGPeaks"), load = TRUE)
