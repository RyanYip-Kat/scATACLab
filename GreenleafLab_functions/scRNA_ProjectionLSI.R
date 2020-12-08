library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(uwot)
library(edgeR)
library(FNN)
library(matrixStats)
library(Rcpp)
library(argparse)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(stringr)
set.seed(1)

functions="/home/ye/Work/R/scATAC/ArchR/GreenleafLab_functions/functions.R"
source(functions)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="./output")


parser$add_argument("--control",
                    type="character",
                    default="0h")

parser$add_argument("--treatment",
                    type="character",
                    default="3h")

parser$add_argument("--ratio",
                    type="double",
                    default=0.75)


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seReference=readRDS(args$control)
seDisease =readRDS(args$treatment)

#Set Clustering Parameters
resolution <- c(0.2,0.8,0.8) #clustering resolution
varGenesToUse <- c(3000,3000,3000) #number of variable genes

print("### OptimizeLSI")
genes.use=intersect(rownames(seReference),rownames(seDisease))
matAll <- cbind(assay(seReference[genes.use,]), assay(seDisease[genes.use,]))
lsiObj <- scRNA_optimizeLSI(matAll, resolution = resolution, varFeatures = varGenesToUse)

#UMAP
print("### Run UMAP")
set.seed(1)
umap <- uwot::umap(
    lsiObj[[length(lsiObj)]]$lsiObj$matSVD[,1:25],
    n_neighbors = 30,
    min_dist = 0.5,
    metric = "euclidean",
    n_threads = 5,
    verbose = TRUE,
    ret_model = FALSE
    )

#Plot Info
cells <- c(rep("reference", ncol(seReference)),rep("disease",ncol(seDisease)))
splitCells <- split(cells,lsiObj[[length(lsiObj)]]$clusters)
df <- data.frame(
    clusters = names(splitCells),
    proportion = unlist(lapply(seq_along(splitCells), function(x) sum(splitCells[[x]]=="disease") / length(splitCells[[x]])))
    )

#Plot UMAP Data Frame
print("### Plot UMAP")
plotDF <- data.frame(umap)
rownames(plotDF) <- c(colnames(seReference), colnames(seDisease))
plotDF$type <- cells
plotDF$clusters <- lsiObj[[length(lsiObj)]]$clusters
plotDF$classification <- 0
#If disease cells are clustered with healthy cluster (proportion > 0.8) we will classify these as healthy-like
plotDF$classification[plotDF$type == "disease" & plotDF$clusters %in% paste0(df$clusters[df[,2] > args$ratio])] <- 1
plotDF$classification[plotDF$type == "disease"] <- plotDF$classification[plotDF$type == "disease"] + 1
plotDF <- plotDF[order(plotDF$classification), ]

#Formal Classification
plotDF$classificationSTR <- "reference"
plotDF$classificationSTR[plotDF$classification==1] <- "healthy-like"
plotDF$classificationSTR[plotDF$classification==2] <- "disease-like"

saveRDS(list(df,plotDF),file.path(args$outdir,"plotDF.rds"))
#Plot PDFs
pdf(file.path(args$outdir,"Classification-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(plotDF, aes(X1,X2,color=classificationSTR)) + 
    geom_point() +
    theme_bw() +
    xlab("UMAP Dimension 1") + 
    ylab("UMAP Dimension 2") +
    scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3"))

dev.off()

################################
#Load Saved UMAP Manifold
#umapManifold <- uwot::load_uwot("data/Supplementary_Data_LSI_Projection/scRNA-Hematopoiesis-UMAP-model.190505.uwot.tar")
print("### Loading umapManifold model")
umap_model=metadata(seReference)$uwot_model
umapManifold=load_uwot(umap_model)

#LSI Projection Matrix
print("### LSI Projection")
lsiGenes <- metadata(seReference)$variableGenes
matProjectLSI <- assay(seDisease[lsiGenes,])

#LSI Project
lsiReference <- metadata(seReference)$optimizeLSI[[length(metadata(seReference)$optimizeLSI)]]$lsiObj
lsiProjection <- projectLSI(matProjectLSI, lsiReference)

#UMAP Projection
#Set Seed Prior to umap_transform (see uwot github)
set.seed(1)
print("### Transform UMAP")
umapProjection <- uwot::umap_transform(as.matrix(lsiProjection), umapManifold, verbose = TRUE)

#Plot Projection
refDF <- data.frame(row.names = colnames(seReference), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference")
proDF <- data.frame(row.names = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = plotDF[colnames(seDisease),]$classificationSTR)
projectionDF <- rbind(refDF, proDF)

pdf(file.path(args$outdir,"Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(projectionDF, aes(X1,-X2,color=Type)) +
    geom_point() +
    theme_bw() +
    xlab("UMAP Dimension 1") +
    ylab("UMAP Dimension 2") +
    scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3"))

dev.off()
saveRDS(umapProjection,file.path(args$outdir,"umapProjection.rds"))
saveRDS(projectionDF,file.path(args$outdir,"projectionDF.rds"))
####################################################
#Differential Analysis Into LSI UMAP
####################################################

print("### Differential Analysis Into LSI UMAP")
#Input Parameters
input_knn <- 25
scaleTo <- 10000
nMax <- 500

#LSI-SVD
svdReference <- as.data.frame(lsiReference$matSVD)
svdDisease <- as.data.frame(as.matrix(lsiProjection))

#Differential Seed
set.seed(1)

#Cells that we are testing of disease
idxDisease <- rownames(plotDF)[plotDF$classificationSTR=="disease-like"]

#If the number of cells exceeds the max downsample to max
if(length(idxDisease) > nMax){
    idxDisease <- sample(idxDisease, nMax)
}

#If the number of cells is greater than 5 continue
stopifnot(length(idxDisease) > 5)

#KNN Nearest Neighbor using FNN
knnDisease <- get.knnx(
    data = svdReference,
    query = svdDisease[idxDisease, ], #Subset by idxDisease
    k = input_knn)

#Determine the minimum KNN where reference cells are less than 1.25x disease cells
i <- 0
uniqueIdx <- unique(as.vector(knnDisease$nn.index))
while(length(uniqueIdx) > 1.25 * length(idxDisease)){
    i <- i + 1
    uniqueIdx <- unique(as.vector(knnDisease$nn.index[,seq_len(input_knn-i)]))
}

#Reference cells for testing
idxReference <- rownames(svdReference)[uniqueIdx]

#If there are more healthy cells downsample healthy cells
#If there are more disease cells downasmple disease cells
if(length(idxReference) > length(idxDisease)){
    idxReference <- sample(idxReference, length(idxDisease))
}else{
    idxDisease <- sample(idxDisease, length(idxReference))
}
message(sprintf("nDisease = %s\nnHealthy = %s", length(idxDisease), length(idxReference)))

#Disease and Reference Matrix
matHealthy <- assay(seReference[,idxReference])
matDisease <- assay(seDisease[,idxDisease])

#Normalize to scaleTo
matNormDisease <- t(t(matDisease)/Matrix::colSums(matDisease)) * scaleTo
matNormHealthy <- t(t(matHealthy)/Matrix::colSums(matHealthy)) * scaleTo

#T-Test Comparisons
print("### Run Comparisons")
dfTT <- sparseMatTTest(matNormDisease, matNormHealthy)
dfTT$feature <- rownames(matNormDisease)
dfTT$log2Mean <- log2(rowMeans(cbind(dfTT$mean1, dfTT$mean2)) + 10^-4)
dfTT$log2FC <- log2((dfTT$mean1 + 10^-4)/(dfTT$mean2 + 10^-4))

plotDiff <- data.frame(row.names=row.names(dfTT),log2Mean=dfTT$log2Mean,log2FC=dfTT$log2FC,FDR=dfTT$fdr)
plotDiff <- plotDiff[complete.cases(plotDiff),]
plotDiff$type <- "not-differential"
plotDiff$type[plotDiff$log2FC > 0.5 & plotDiff$FDR < 0.01] <- "up-regulated"
plotDiff$type[plotDiff$log2FC < -0.5 & plotDiff$FDR < 0.01] <- "do-regulated"

pdf(file.path(args$outdir,"Differential-MA-Plot.pdf"), width = 8, height = 6, useDingbats = FALSE)
ggplot(plotDiff, aes(log2Mean,log2FC,color=type)) +
    geom_point(size=0.5) +
    theme_bw() +
    xlab("log2 Mean") +
    ylab("log2 Fold Change") +
    scale_color_manual(values=c("not-differential"="lightgrey", "do-regulated"="dodgerblue3", "up-regulated"="firebrick3"))

dev.off()

#Save Output
readr::write_tsv(dfTT, file.path(args$outdir,"Differential-Results.tsv"))

