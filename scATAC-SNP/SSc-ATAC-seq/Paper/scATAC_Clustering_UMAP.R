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

####################################################
#For Clustering Analysis Start Here
####################################################
nPCs1 <- 1:50 #Number of PCs in first analysis using all peaks
nTop <- 50000 #Choose a higher number of variable peaks across clusters (25-50k) to mitigate batch effects
nPCs2 <- 1:50 #Number of PCs in second analysis using variable peaks across clusters
resolution <- 1.5 #Clustering resolution for Seurat SNN

#RUN LSI 1
message("Running LSI 1...")
mat <- assay(se)
lsi1 <- calcLSI(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL)
clust1 <- seuratSNN(lsi1[[1]], dims.use = nPCs1, resolution = resolution)

#Make Pseudo Bulk Library
message("Making PseudoBulk...")
mat <- mat[,rownames(lsi1[[1]]), drop = FALSE] #sometimes cells are filtered
mat@x[mat@x > 0] <- 1 #binarize
clusterSums <- groupSums(mat = mat, groups = clust1, sparse = TRUE) #Group Sums
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), nTop) #Top variable peaks

#RUN LSI 2
message("Running LSI 2...")
lsi2 <- calcLSI(mat[varPeaks,,drop=FALSE], nComponents = 50, binarize = TRUE, nFeatures = NULL)
clust2 <- seuratSNN(lsi2[[1]], dims.use = nPCs2, resolution = resolution)

#Append Summarized Experiment
se <- se[,rownames(lsi2[[1]])]
colData(se)$Clusters <- clust2
metadata(se)$LSI <- lsi2
metadata(se)$LSIPeaks <- varPeaks
metadata(se)$matSVD <- lsi2$matSVD


#Set Seed and perform UMAP on LSI-SVD Matrix
matSVD <- metadata(se)$matSVD
clusters <- colData(se)$Clusters

set.seed(1)
message("INFO : Run UMAP")
uwotUmap <- uwot::umap(
    matSVD[,1:50],
    n_neighbors = 55,
    min_dist = 0.45,
    metric = "euclidean",
    n_threads = 1,
    verbose = TRUE,
    ret_nn = TRUE,
    ret_model = TRUE
    )

pdf(file.path(outDir,"Plot_UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
df <- data.frame(
    x = uwotUmap[[1]][,1],
    y = uwotUmap[[1]][,2],
    color = clusters
    )
ggplot(df,aes(x,y,color=color)) +
    geom_point() +
    theme_bw() +
    #scale_color_manual(values=metadata(se)$colorMap$Clusters) +
    xlab("UMAP Dimension 1") +
    ylab("UMAP Dimension 2")
dev.off()

colData(se)$UMAP1 <- uwotUmap[[1]][,1]
colData(se)$UMAP2 <- uwotUmap[[1]][,2]

metadata(se)$UMAP_Params <- list(NN = 55, MD = 0.45, PCs = 1:50, VarPeaks = 50000, Res = "1.5")
message("INFO : Save")
saveRDS(se,file.path(outDir,"scATAC-Summarized-Experiment-UMAP.rds"))
