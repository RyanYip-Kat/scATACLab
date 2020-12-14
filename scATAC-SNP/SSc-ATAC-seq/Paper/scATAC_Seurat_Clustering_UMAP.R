library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
library(Seurat)
library(stringr)
library(argparse)
library(ggplot2)

source("/home/ye/Work/BioAligment/SNP/Shi/SSc-ATAC-seq/Paper/functions.R")
################# functions


####################################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--se",
                    type="character",
                    default=NULL,
                    help="Summarized-Experiment object")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
outDir=args$outdir
makedir(outDir)

message("INFO : Loading dataset")
se=readRDS(args$se)


nTop <- 25000
nPCs1 <- 1:50
nPCs2 <- 1:50

message("Making Seurat LSI Object...")
obj <- seuratLSI(assay(se), nComponents = max(nPCs1), nFeatures = NULL)
stopifnot(identical(rownames(obj@meta.data), colnames(se)))
obj@meta.data <- as.data.frame(cbind(obj@meta.data, colData(se)))

message("Adding Graph Clusters...")
obj <- FindNeighbors(obj,reduction="pca",dims =nPCs1)
obj <- FindClusters(object = obj)

#Make Pseudo Bulk Library
mat <- assay(se)
mat@x[mat@x > 0] <- 1
clusterSums <- groupSums(mat = mat, groups = paste0("C",obj$seurat_clusters), sparse = TRUE)
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3)
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), nTop)

#Re-run Seurat LSI
message("Making Seurat LSI Object...")
obj2 <- seuratLSI(assay(se)[varPeaks,], nComponents = max(nPCs2), nFeatures = NULL)
stopifnot(identical(rownames(obj2@meta.data), colnames(se)))
obj2@meta.data <- as.data.frame(cbind(obj2@meta.data, colData(se)))

message("Adding Graph Clusters...")
obj2 <- FindNeighbors(obj2,reduction="pca",dims =nPCs2)
obj2 <- FindClusters(object = obj2)

#Plot uMAP
message("Running UMAP")
obj2 <- RunUMAP(object = obj2,reduction = "pca", dims= nPCs2)
plotUMAP <- data.frame(Embeddings(obj2,reduction="umap"), obj2@meta.data)
colnames(plotUMAP) <- c("x","y",colnames(plotUMAP)[3:ncol(plotUMAP)])
colData(se)$Clusters <- paste0("Cluster",as.integer(plotUMAP[,"seurat_clusters"]) + 1)
colData(se)$UMAP1 <- plotUMAP$x
colData(se)$UMAP2 <- plotUMAP$y

message("INFO : Save")
saveRDS(se, file.path(outDir,"scATAC-Summarized-Experiment.rds"))

pdf(file.path(outDir,"LSI-Clustering-Peaks.pdf"))
ggplot(plotUMAP, aes(x=x,y=y,color=seurat_clusters)) + geom_point(size = 0.5) + 
  theme_bw() + xlab("UMAP1") + ylab("UMAP2")
dev.off()

