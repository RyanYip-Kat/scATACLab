library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(rtracklayer)
library(readr)
library(GenomicRanges)
library(magrittr)
library(stringr)
library(edgeR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(argparse)

makedir=function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}



########################################   callpeak paramenters
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--name",
                    type="character",
                    default="CEMT",
                    help="name add to gr mcols")

parser$add_argument("--meta",
                    type="character",
                    default=NULL,
                    help="*.bw file in meta csv file")


args <- parser$parse_args()

########################################
message("Write Paper bw file into BED files")
DF=read.csv(args$meta,header=FALSE,sep=",",stringsAsFactors=FALSE)
bwFiles=as.character(DF$V2)
Names=as.character(DF$V1)


outPath=args$outdir
makedir(outPath)


bedList=list()
for(i in 1:length(bwFiles)){
   name=Names[i]
   bw=bwFiles[i]
   cat(sprintf("INFO : [ %d of %d ] --- [ %s ] \n",i,length(bwFiles),bw))
   gr=import(bw)
   mcols(gr)$Sample=name
   bed=as.data.frame(gr)
   out=bed[,c("seqnames","start","end","Sample")]
   outfile=file.path(outPath,paste0(name,".bed"))
   write.table(bed[,c("seqnames","start","end")],outfile,sep="\t",row.names=F,col.names=F)
   bedList[[name]]=out
}

AllBed=do.call(rbind,bedList)
gr=makeGRangesFromDataFrame(AllBed)
mcols(gr)$Disease=args$name

saveRDS(gr,file.path(outPath,"GR.rds"))
