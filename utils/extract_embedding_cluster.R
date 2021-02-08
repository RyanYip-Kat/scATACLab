library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--archr",
                    type="character",
		    default=NULL,
                    help="the path of project saved")

parser$add_argument("--column",
                    type="character",
                    default="Clusters",
		    help="the column in cellcol to be exported")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

parser$add_argument("--loupe_cell",
                    action="store_true",
                    default=FALSE)

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading data")
proj<-readRDS(args$archr)


print("# Get cellcolData from ArchR Project")
cells=getCellNames(proj)
if(args$loupe_cell){
	cells=proj$Barcode
}
mat=getEmbedding(proj,"TSNE")
colnames(mat)=c("tSNE_1","tSNE_2")
mat$Barcode=cells
mat=mat[,c(3,1,2)]
write.table(mat,file.path(args$outdir,"tsne_projection.csv"),sep=",",quote=F,row.names = F)

mat=getEmbedding(proj,"UMAP")
colnames(mat)=c("UMAP_1","UMAP_2")
mat$Barcode=cells
mat=mat[,c(3,1,2)]
write.table(mat,file.path(args$outdir,"umap_projection.csv"),sep=",",quote=F,row.names = F)

meta=as.data.frame(getCellColData(proj,select=args$column))
meta$barcode=cells
meta=meta[,c(2,1)]
write.table(meta,file.path(args$outdir,"cluster.csv"),sep=",",quote=F,row.names = F)


