library(Seurat)
library(argparse)
library(stringr)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--reduction",
                    type="character",
                    default=NULL,help="tsne,umap")

parser$add_argument("--seurat",
                    type="character",
                    default="")
args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

object=readRDS(args$seurat)
print("#### Add reduction")
mat=read.csv(args$reduction,stringsAsFactors=F)
rownames(mat)=mat$Barcode
mat=mat[colnames(object),]
mat=mat[,-1]
pcs=colnames(mat)
mat=as.matrix(mat)
key=str_split(colnames(mat)[1],"_")[[1]][1]
r=str_to_lower(key)
object[[r]]<-CreateDimReducObject(embeddings =mat,
                                  key =paste0(key,"_"),
                                  assay = DefaultAssay(object))

saveRDS(object,file.path(args$outdir,"seurat.rds"))
