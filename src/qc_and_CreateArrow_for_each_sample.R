library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--inputFiles",
                    type="character",
                    default="")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--genome",
                    type="character",
                    default="mm10",
		    choices=c("hg38","hg19","mm9","mm10"),
		    help="which genome to be used")

parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")


parser$add_argument("--rm_doublet",
		    action="store_true",
                    default=FALSE,
                    help="whether to remove doublet")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


inputFiles<-read.table(args$inputFiles,sep=",",stringsAsFactors=FALSE)
files=as.character(inputFiles$V1)
SampleNames=as.character(inputFiles$V2)
Tss=as.numeric(inputFiles$V3)
Frags=as.numeric(inputFiles$V4)

print(paste0("Setting default genome  with :",args$genome))
addArchRGenome(args$genome)

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 


###################
max_frags=50000
max_tss=15

#################
print("# Create ArrowFiles")
ArrowFiles=c()
for(i in 1:length(files)){
	arrow=createArrowFiles(
			       inputFiles = files[i],
                               sampleNames =SampleNames[i],
                               filterTSS = Tss[i], #Dont set this too high because you can always increase later
                               filterFrags = Frags[i],
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
			       minFrags = 500,
			       maxFrags = max_frags,
			      )
	ArrowFiles=c(ArrowFiles,arrow)
	print(arrow)
}
print(ArrowFiles)

projHeme <- ArchRProject(
  ArrowFiles = ArrowFiles,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

if(args$rm_doublet){
        projHeme <- addDoubletScores(
                    input = projHeme,
                    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
                    LSIMethod = 1)

        projHeme <- filterDoublets(projHeme)
}

print("### Save ArchR Object")
saveArchRProject(ArchRProj = projHeme, outputDirectory =file.path(args$outdir,"Save-ProjHeme-Create"), load = TRUE)
