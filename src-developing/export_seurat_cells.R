library(argparse)
library(stringr)
library(Seurat)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--seurat",
                    type="character",
                    default=NULL)

parser$add_argument("--rename",
                    action="store_true",
                    default=FALSE)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seurat=readRDS(args$seurat)

counts=GetAssayData(seurat,"counts")
if(args$rename){
	cellNames=unlist(lapply(colnames(counts),function(cell)return(str_split(cell,"-")[[1]][1])))
	cellNames=paste(cellNames,"1",sep="-")
	colnames(counts)=cellNames
}

cells=data.frame("Barcode"=colnames(counts))
genes=data.frame("Gene"=str_to_title(rownames(counts)))

write.table(cells,file.path(args$outdir,"barcode.csv"),quote=FALSE,sep=",",row.names=FALSE)
write.table(genes,file.path(args$outdir,"genes.csv"),quote=FALSE,sep=",",row.names=FALSE)

