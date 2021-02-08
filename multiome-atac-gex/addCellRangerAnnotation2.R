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
		    nargs="+",
                    type="character",
                    default=NULL,
		    help="*atac_peak_annotation.tsv")


parser$add_argument("--annoName",
		    nargs="+",
		    type="character",
		    default=NULL,
		    help="name for each annotation")


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

AnnoFiles=args$annotation
names(AnnoFiles)=args$annoName
message("INFO : Loading Annotation ...")
keepChr=getSeqnames(proj)
AnnoList=lapply(seq_along(AnnoFiles),function(i){
			AnnoDF=read.table(AnnoFiles[i],sep="\t",header=TRUE)
			AnnoDF=subset(AnnoDF,str_detect(peak,paste(keepChr,collapse="|")))
			cat(sprintf("INFO : make [ %s ] Annotation Granges ...\n",names(AnnoFiles)[i]))
			DFList=lapply(AnnoDF$peak,function(peak){
					      pl=str_split(peak,"_|-")[[1]]
					      x=data.frame("chr"=pl[1],"start"=pl[2],"end"=pl[3])
					      return(x)})
			peakDF=do.call(rbind,DFList)
			AnnoDF=cbind(peakDF,AnnoDF)
			AnnoGR=makeGRangesFromDataFrame(AnnoDF,keep.extra.columns=TRUE)
			return(AnnoGR)
		    })

names(AnnoList)=names(AnnoFiles)
AnnoGRList=as(AnnoList,"GRangesList")
saveRDS(AnnoGRList,file.path(args$outdir,"GRangesList.rds"))
message("INFO : addPeakAnnotations ...")
proj=addPeakAnnotations(proj,regions=AnnoGRList,name="CellRanger")

message("INFO : Save ..")
saveArchRProject(ArchRProj = proj,overwrite=TRUE,outputDirectory =file.path(args$outdir,"Save-ProjHeme-CellRangerAnno"), load = TRUE)
