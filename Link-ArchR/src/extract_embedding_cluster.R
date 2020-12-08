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

meta=as.data.frame(getCellColData(proj))
write.table(meta,file.path(args$outdir,"metadata.csv"),sep=",",quote=F,row.names = F)


