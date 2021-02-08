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
		    default=NULL,
                    help="the path of project saved")

parser$add_argument("--columns",
		    nargs="+",
                    type="character",
                    default=NULL,
		    help="the columns in cellcol to be exported")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")


args <- parser$parse_args()

options(stringsAsFactors=FALSE)
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading data")
proj=loadArchRProject(args$project)


print("# Get cellcolData from ArchR Project")
cells=proj$Barcode
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

metadata=as.data.frame(getCellColData(proj))
for(column in args$columns){
	cat(sprintf("### Export Column : %s \n",column))
	meta=data.frame("barcode"=cells,column=metadata[[column]])
	colnames(meta)=c("barcode",column)
        write.table(meta,file.path(args$outdir,paste0(column,"_cluster.csv")),sep=",",quote=F,row.names = F)
}



