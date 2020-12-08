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


parser$add_argument("--column",
                    type="character",
                    default="Clusters",
                    help="the column in metadata for cross table")


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

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)

meta=as.data.frame(getCellColData(projHeme))
if(!args$column%in%colnames(meta)){
	stop("the input column not in project's meta,please run update_metadata first!!!")
}

stopifnot("Sample"%in%colnames(meta))
m=table(meta[["Sample"]],meta[[args$column]])
write.table(as.matrix(m),file.path(args$outdir,paste0("Sample_",args$column,"_number.csv")),sep=",",quote=F)

stopifnot("Status"%in%colnames(meta))
m=table(meta[["Status"]],meta[[args$column]])
write.table(as.matrix(m),file.path(args$outdir,paste0("Status_",args$column,"_number.csv")),sep=",",quote=F)
