library(Seurat)
library(Signac)
library(argparse)
library(stringr)
library(future)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")

parser$add_argument("--groupby",
                    type="character",
                    default="seurat_clusters",
		    help="group in metadata used for callpeaking")


args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

findMacs2=function(){
	macs2_path=system("which macs2",intern=TRUE)
	return(macs2_path)
}

############################## Configure
plan("multiprocess", workers = 16)
outDir=dirname(args$seurat)
callpeakDir=file.path(outDir,"CallPeak")
Group=args$groupby

makedir(callpeakDir)
############################## Loading dataset
message("INFO : Loading dataset")
object=readRDS(args$seurat)


############################## callpeaks
cat(sprintf("INFO : callpeak with [ %s ]\n",Group))
peaks <- CallPeaks(
  object = object,
  group.by = Group,
  macs2.path = findMacs2(),
  outdir=callpeakDir,
  cleanup=FALSE,
  effective.genome.size=2.7e+09      # human
)

message("INFO : Save Object")
saveRDS(peaks,file.path(outDir,"callpeaks.rds"))
