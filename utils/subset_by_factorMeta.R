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
		    help="the column in cellcol to be exported")

parser$add_argument("--subset",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset of column to be extracted")

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
projHeme=loadArchRProject(args$project)

print("# Get cellcolData from ArchR Project")
metadata=as.data.frame(getCellColData(projHeme))
cells=getCellNames(projHeme)


if(is.null(args$column) | is.null(args$subset)){
	stop("column and subset for subset archr can not be null!!!")
}	
stopifnot(args$column%in%colnames(metadata))
target=metadata[[args$column]]
idxPass= which(target%in%args$subset)

stopifnot(args$subset%in%unique(target))

cellPass=cells[idxPass]
subsetArchRProject(projHeme,cells=cellPass,outputDirectory=file.path(args$outdir,"ArchRSubset"))

	
