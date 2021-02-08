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
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}




print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
proj<-loadArchRProject(args$project)

arrowFiles=getArrowFiles(proj)
scores <- addDoubletScores(
                    input = arrowFiles,
                    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
                    LSIMethod = 1)


scoreDF=list()
for(name in names(scores)){
	d=data.frame("DoubletScore"=scores[[name]][["doubletScore"]],"DoubletEnrichment"=scores[[name]][["doubletEnrich"]])
	#rownames(d)=names(scores[[name]][["doubletScore"]])
	scoreDF[[name]]=d
}

DF=Reduce(rbind,scoreDF)
DF=DF[getCellNames(proj),]

write.table(DF,file.path(args$outdir,"DoubletScores.csv"),sep=",",quote=FALSE)
for(name in colnames(DF)){
        cat(sprintf("INFO : %s \n",name))
        proj=addCellColData(proj,data=DF[,name],name=name,cells=rownames(DF))
    }

proj=filterDoublets(proj)
saveRDS(proj,file.path(getOutputDirectory(proj),"Save-ArchR-Project.rds"))


