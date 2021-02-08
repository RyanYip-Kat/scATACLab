library(argparse)
library(stringr)
#library(Seurat)
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

parser$add_argument("--harmony",
                    action="store_true",
                    default=FALSE,
                    help="wether to runn harmony")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting default genome  with :",args$genome))
addArchRGenome(args$genome)

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)

print("# Run LSI ")
projHeme <- addIterativeLSI(
    ArchRProj = projHeme,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 10,
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2),
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 25000,
    dimsToUse = 1:30,
    force=TRUE
)

if(args$harmony){
	print("# Run harmony ")
	projHeme <- addHarmony(
			       ArchRProj = projHeme,
			       reducedDims = "IterativeLSI",
                               name = "Harmony",
                               groupBy = "Sample",
			       force=TRUE)
}


use_rep=ifelse(args$harmony,"Harmony","IterativeLSI")

print("# Run tSNE ")
projHeme <- addTSNE(
    ArchRProj = projHeme,
    reducedDims =use_rep,
    name="TSNE",
    method="RTSNE",
    force=TRUE
)

print("# Run UMAP ")
projHeme <- addUMAP(
    ArchRProj = projHeme,
    reducedDims =use_rep,
    name = "UMAP",
    force=TRUE
)

print("# Run Clusters ")
projHeme <- addClusters(
    input = projHeme,
    reducedDims = use_rep,
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force=TRUE
)
print(table(projHeme$Clusters))

projHeme <- addImputeWeights(projHeme)
cM <- confusionMatrix(paste0(projHeme$Clusters), paste0(projHeme$Sample))
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
pdf(file.path(args$outdir,"pheatmap.pdf"),width=12,height=16)
pheatmap(
    mat = as.matrix(cM),
    color = viridis(100),
    border_color = "black"
)
dev.off()

print("# Save ...")
saveArchRProject(ArchRProj = projHeme, outputDirectory =file.path(args$outdir,"Save-ProjHeme"), load = TRUE)
