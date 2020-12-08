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

parser$add_argument("--metadata",
                    type="character",
                    default=NULL,
		    help="the metadata to be add on project")

parser$add_argument("--name",
                    type="character",
                    default="new_cluster",
                    help="the metadata name to be add on project")


#parser$add_argument("--outdir",
#                    type="character",
#                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
#if(!dir.exists(args$outdir)){
#        dir.create(args$outdir,recursive=TRUE)
#}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)

meta=as.data.frame(getCellColData(projHeme))
meta$Barcode=getCellNames(projHeme)
#if(!"Barcode"%in%colnames(meta)){
#	stop("Barcode not in project's meta,please run update_metadata first!!!")
#}

Data=read.csv(args$metadata,stringsAsFactors=FALSE,header=TRUE,row.names=1)  # barcode,cluster(GGAACCCCAGGTGGTA-8,"C1")
stopifnot(ncol(Data)==1)
colnames(Data)=c("cluster")
Data=Data[meta$Barcode,,drop=FALSE]

print("# Add  data on CellColData")
if(args$name%in%colnames(meta)){
	name=paste0("new_",args$name)
}else{
	name=args$name
}
print(paste0("## Add :",name," on CellColData"))
print(head(Data$cluster))
projHeme=addCellColData(projHeme,data=Data$cluster,name=name,cells=rownames(meta))
print("# Save ...")
path=getOutputDirectory(projHeme)
saveArchRProject(ArchRProj = projHeme, outputDirectory =path, load = TRUE)


