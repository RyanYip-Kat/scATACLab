library(Seurat)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
		    type="character",
		    default="")

parser$add_argument("--column",
                    default="Clusters",
                    type="character")

parser$add_argument("--name",
                    default=NULL,
                    type="character")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)

mat<-seurat@meta.data[,args$column,drop=FALSE]
mat$barcode=rownames(mat)
mat<-mat[,c(2,1)]
if(!is.null(args$name)){
	mat=data.frame("barcode"=mat$barcode,"cluster"=args$name)
	colnames(mat)=c("barcode",args$column)
}
write.table(mat,file.path(args$outdir,paste0(args$column,"_clusters.csv")),sep=",",quote=FALSE,row.names=FALSE)


