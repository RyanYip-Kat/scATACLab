library(Seurat)
library(argparse)
library(Signac)
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
                    default=NULL,
                    type="character")

#parser$add_argument("--name",
#                    default="BC",
#                    type="character")

parser$add_argument("--reduction",
		    action="store_true",
		    default=FALSE)
args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat<-readRDS(args$seurat)
if(args$reduction){
	#mat<-Embeddings(seurat,"tsne")
        #mat<-as.data.frame(mat)
        #mat$Barcode<-rownames(mat)
        #mat<-mat[,c(3,1,2)]
        #write.table(mat,file.path(args$outdir,"tsne_projection.csv"),sep=",",quote=FALSE,row.names=FALSE)

	mat<-Embeddings(seurat,"umap")
        mat<-as.data.frame(mat)
        mat$Barcode<-rownames(mat)
        mat<-mat[,c(3,1,2)]
        write.table(mat,file.path(args$outdir,"umap_projection.csv"),sep=",",quote=FALSE,row.names=FALSE)
}

if(!is.null(args$column)){
	mat<-seurat@meta.data[,args$column,drop=FALSE]
        mat$barcode=rownames(mat)
        mat<-mat[,c(2,1)]
        write.table(mat,file.path(args$outdir,paste0(args$column,"_clusters.csv")),sep=",",quote=FALSE,row.names=FALSE)
}

#mat=data.frame("barcode"=colnames(seurat),"name"=rep(args$name,ncol(seurat)))
#colnames(mat)=c("barcode","name")
#write.table(mat,file.path(args$outdir,"name_clusters.csv"),sep=",",quote=FALSE,row.names=FALSE)



