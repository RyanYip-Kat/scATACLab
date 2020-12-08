library(argparse)
library(stringr)
library(future)
library(Signac)
library(Seurat)
library(patchwork)
library(ggplot2)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")

parser$add_argument("--peak",
                    type="character",
                    default=NULL,
                    help="call peak result from (CallPeaks)")


parser$add_argument("--region",
		    nargs="+",
                    type="character",
                    default="CD8A",
		    help="genes name to be plotted")

parser$add_argument("--groupby",
                    type="character",
                    default="seurat_clusters",
                    help="group use for CoveragePlot")



parser$add_argument("--outdir",
                    type="character",
                    default="./Results")

args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

################################
outDir=file.path(args$outdir,"CoveragePlot")
makedir(outDir)


################################
message("INFO : Loading  dataset ...")
seurat=readRDS(args$seurat)
peaks=readRDS(args$peak)


################################ plot regions
regions=args$region
for(i in seq_along(regions)){
	region=regions[i]
	cat(sprintf("INFO : CoveragePlot [ %d of %d ]  --- [ %s ]\n",i,length(regions),region))
	tryCatch({ p=CoveragePlot(
		       object = seurat,
		       region = region,
		       ranges = peaks,
		       ranges.title = "MACS2")
	           fileName=file.path(outDir,paste0(region,".pdf"))
	           ggsave(fileName,plot=p,width=16,height=10)},
		   error=function(e){cat(sprintf("INFO : Invalid region : [ %s ]\n",region))},
		   finnally={cat(sprintf("INFO : CoveragePlot [ %d of %d ]  --- [ %s ] sucessfuly!\n",i,length(regions),region))}
	)
}



