library(argparse)
library(stringr)
library(Seurat)
library(Signac)
library(ArchR)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(viridis)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--genome",
                    type="character",
                    default="hg38",
		    choices=c("hg38","hg19","mm9","mm10"),
		    help="which genome to be used")

parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("### Setting default genome  with :",args$genome))
addArchRGenome(args$genome)

print(paste0("### Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("### Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)

matrixs=c("PeakMatrix","GeneScoreMatrix")
avaliable_matrixs=getAvailableMatrices(projHeme)
print(avaliable_matrixs)

avaliable_matrixs=matrixs[matrixs%in%avaliable_matrixs]
metadata=as.data.frame(getCellColData(projHeme))
stopifnot("Barcode"%in%colnames(metadata))

data_list<-list()
for(mat in avaliable_matrixs){
	print(paste0("### Get Matrix from ",mat," slot"))
	binary=ifelse(mat=="TileMatrix",TRUE,FALSE)
	se=getMatrixFromProject(projHeme,useMatrix=mat,binarize=binary)
	features=str_to_upper(getFeatures(projHeme,useMatrix=mat))
	counts=assay(se)
	colnames(counts)=colData(se)$Barcode
	rownames(counts)=features
	data_list[[mat]]=counts
}

data_list[["metadata"]]=metadata
print("### Save Data List")
saveRDS(data_list,file.path(args$outdir,"ArchR_DataList.rds"))
