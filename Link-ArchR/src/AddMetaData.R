library(argparse)
library(stringr)
library(Seurat)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default="")

parser$add_argument("--metadata",
                    type="character",
                    default=NULL)

parser$add_argument("--name",
                    type="character",
                    default=NULL)
args <- parser$parse_args()

print("### Subset")
DATA=read.csv(args$metadata,header=TRUE,stringsAsFactors=F)
colnames(DATA)=c("barcode","celltype")
rownames(DATA)<-DATA$barcode


DATA<-subset(DATA,select=celltype)
seurat<-readRDS(args$seurat)

if(nrow(DATA)!=ncol(seurat)){
	stop("metadata dim must be equal to seurat")
}

print("Update Metadata")
if(is.null(args$name)){
	name="celltype"
}else{
	name=args$name
}
seurat<-AddMetaData(seurat,metadata=DATA,col.name=name)
print("Save")
saveRDS(seurat,args$seurat)
