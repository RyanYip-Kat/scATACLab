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


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

parser$add_argument("--batch_correct",
		    action="store_true")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



#print(paste0("Setting default genome  with :",args$genome))
#addArchRGenome(args$genome)

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
proj<-loadArchRProject(args$project)
#if(is.null(getImputeWeights(proj))){
#        proj=addImputeWeights(proj)
#}

message("INFO : Run LSI-ATAC ...")
proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(
      resolution = 0.8,
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    name = "LSI_ATAC",
    force=TRUE
)


message("INFO : Run LSI-RNA ...")
proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(
      resolution = 0.8,
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix",
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
)

message("INFO : Combined Dims ...")
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

if(is.null(getImputeWeights(proj))){
        message("INFO : imputeWeigth ...")
        proj=addImputeWeights(proj,reducedDims="LSI_Combined")
}


RNA_rep="LSI_RNA"
ATAC_rep="LSI_ATAC"
Combined_rep="LSI_Combined"
if(args$batch_correct){
	message("INFO :  Run harmony : LSI_ATAC ")
        proj <- addHarmony(
			   ArchRProj = proj,
                           reducedDims = "LSI_ATAC",
                           name = "Harmony_ATAC",
                           groupBy = "Sample",
                           force = TRUE)

        message("INFO :  Run harmony : LSI_RNA ")
        proj <- addHarmony(
			   ArchRProj = proj,
                           reducedDims = "LSI_RNA",
                           name = "Harmony_RNA",
                           groupBy = "Sample",
                           force = TRUE)
	
	message("INFO :  Run harmony : LSI_Combined ")
        proj <- addHarmony(
			   ArchRProj = proj,
			   reducedDims = "LSI_Combined",
                           name = "Harmony_Combined",
			   groupBy = "Sample",
			   force = TRUE)
	RNA_rep="Harmony_RNA"
	ATAC_rep="Harmony_ATAC"
	Combined_rep="Harmony_Combined"
}


message("INFO : Run ATAC UMAPs ...")
proj <- addUMAP(proj, reducedDims = ATAC_rep, name = "UMAP_ATAC",  force = TRUE)

message("INFO : Run RNA UMAPs ...")
proj <- addUMAP(proj, reducedDims = RNA_rep, name = "UMAP_RNA",  force = TRUE)

message("INFO : Run Combined UMAPs ...")
proj <- addUMAP(proj, reducedDims = Combined_rep, name = "UMAP_Combined",  force = TRUE)

message("INFO : Add Clusters ...")
proj <- addClusters(proj, reducedDims = Combined_rep, name = "Clusters", resolution = 0.8, force = TRUE)

message("INFO : Save ArchR Object ...")
saveArchRProject(ArchRProj = proj, outputDirectory =file.path(args$outdir,"Save-ProjHeme-Clusters"), load = TRUE)
saveRDS(proj,file.path(args$outdir,"Save-ProjHeme-Clusters/Save-ArchR-Project.rds"))

message("INFO : Plot demo ...")
p1 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)
plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)
message("INFO : Done!")
