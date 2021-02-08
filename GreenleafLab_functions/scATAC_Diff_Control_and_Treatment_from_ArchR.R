library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(uwot)
library(edgeR)
library(FNN)
library(matrixStats)
library(Rcpp)
library(argparse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
set.seed(1)

functions="/home/ye/Work/R/scATAC/ArchR/GreenleafLab_functions/functions.R"
source(functions)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="./output")


parser$add_argument("--combinedSE",
                    type="character",
                    default=NULL,
                    help="the control se object")

parser$add_argument("--control_se",
                    type="character",
                    default=NULL)

parser$add_argument("--treatment_se",
                    type="character",
                    default=NULL)

parser$add_argument("--control_umap",
                    type="character",
                    default=NULL)

parser$add_argument("--promoterPeaks",
                    type="character",
                    default=NULL)

parser$add_argument("--plotDF",
                    type="character",
                    default=NULL)


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
combinedSE=readRDS(args$combinedSE)
plotDF=readRDS(args$plotDF)

seControl <- readRDS(args$control_se)
seDisease=readRDS(args$treatment_se)
#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot(args$control_umap)

#LSI Projection Matrix
#lsiPeaks <- metadata(se)$variablePeaks
#matProjectLSI <- assay(seDisease[lsiPeaks,])

#LSI Project
lsiReference <- metadata(seControl)$IterativeLSI
#lsiProjection <- projectLSI(matProjectLSI, lsiReference)
lsiDisease=metadata(seDisease)$IterativeLSI

#UMAP Projection
#Set Seed Prior to umap_transform (see uwot github)
set.seed(1)
umapProjection <- uwot::umap_transform(as.matrix(lsiDisease$matSVD), umapManifold, verbose = TRUE)

#Plot Projection
refDF <- data.frame(row.names = colnames(seControl), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference")
proDF <- data.frame(row.names = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = plotDF[colnames(seDisease),]$classificationSTR)
projectionDF <- rbind(refDF, proDF)

pdf(file.path(args$outdir,"Projection-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)
ggplot(projectionDF, aes(X1,X2,color=Type)) +
    geom_point() +
    theme_bw() +
    xlab("UMAP Dimension 1") +
    ylab("UMAP Dimension 2") +
    scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3"))

dev.off()

####################################################
#Differential Analysis Into LSI UMAP
####################################################

#Previous MPAL and Reference Summarized Experiment
#Contains Peaks for MPALs and Reference Cell Union Set
promoterPeaks <- readRDS(args$promoterPeaks)

#Input Parameters
input_knn <- 25
scaleTo <- 10000
nMax <- 500

#LSI-SVD
svdReference <- as.data.frame(lsiReference$matSVD)
svdDisease <- as.data.frame(as.matrix(lsiDisease$matSVD))

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
matHealthy <- assay(combinedSE[,idxReference])
matDisease <- assay(combinedSE[,idxDisease])

#Normalize to scaleTo
matNormDisease <- t(t(matDisease)/Matrix::colSums(matDisease[promoterPeaks,])) * 5000
matNormHealthy <- t(t(matHealthy)/Matrix::colSums(matHealthy[promoterPeaks,])) * 5000

#T-Test Comparisons
dfTT <- sparseMatTTest(matNormDisease, matNormHealthy)
dfTT$feature <- rownames(matNormDisease)
dfTT$log2Mean <- log2(rowMeans(cbind(dfTT$mean1, dfTT$mean2)) + 10^-4)
dfTT$log2FC <- log2((dfTT$mean1 + 10^-4)/(dfTT$mean2 + 10^-4))

plotDiff <- data.frame(row.names=row.names(dfTT),log2Mean=dfTT$log2Mean,log2FC=dfTT$log2FC,FDR=dfTT$fdr)
plotDiff <- plotDiff[complete.cases(plotDiff),]
plotDiff$type <- "not-differential"
plotDiff$type[plotDiff$log2FC >= 0.5 & plotDiff$FDR <= 0.05] <- "up-regulated"
plotDiff$type[plotDiff$log2FC <= -0.5 & plotDiff$FDR <= 0.05] <- "do-regulated"


pdf(file.path(args$outdir,"Differential-MA-Plot.pdf"), width = 8, height = 6, useDingbats = FALSE)
ggplot(plotDiff, aes(log2Mean,log2FC,color=type)) +
    #ggrastr::geom_point_rast(size=0.5) +
    geom_point()+
    theme_bw() +
    xlab("log2 Mean") +
    ylab("log2 Fold Change") +
    scale_color_manual(values=c("not-differential"="lightgrey", "do-regulated"="dodgerblue3", "up-regulated"="firebrick3"))

dev.off()

#Save Output
readr::write_tsv(dfTT, file.path(args$outdir,"Differential-Results.tsv"))

