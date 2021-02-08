library(argparse)
library(stringr)
#library(Seurat)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--inputFiles",
                    type="character",
                    default="")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


inputFiles<-read.table(args$inputFiles,sep=",",stringsAsFactors=FALSE)
files=as.character(inputFiles$V1)
filenames=as.character(inputFiles$V2)
barcode_csv=as.character(inputFiles$V3)

print("# Get valid barcode")
lapply(1:length(barcode_csv),function(i){
	       csv=barcode_csv[i]
	       d=read.csv(csv,stringsAsFactors=FALSE,header=TRUE)
	       d=subset(d,select=Barcode)
	       colnames(d)="barcode"
	       d$cell_id=1:nrow(d)
	       write.table(d,file.path(args$outdir,paste0(filenames[i],"_validBarcode.csv")),sep=",",quote=F,row.names=F)})

