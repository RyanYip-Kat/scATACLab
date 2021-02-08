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


parser$add_argument("--barcode",
                    type="character",
                    default=NULL,
		    help="the barcode csv file")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=4,
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

print("# Loading barcode file")
Data=read.csv(args$barcode,stringsAsFactors=FALSE)
colnames(Data)=c("cells","celltype")
cells=as.character(Data$cells)

target=metadata[["Barcode"]]  #consistent with loupe file 
idxPass= which(target%in%cells)


cellPass=cells[idxPass]
subsetArchRProject(projHeme,cells=cellPass,outputDirectory=file.path(args$outdir,"ArchRSubset"))

	
