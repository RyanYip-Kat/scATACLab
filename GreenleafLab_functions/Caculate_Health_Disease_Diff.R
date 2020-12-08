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
                    default="ArchR_result")


parser$add_argument("--treanment",
                    type="character",
                    default=NULL,
                    help="the treanment se object")

parser$add_argument("--control",
                    type="character",
                    default=NULL,
		    help="the control se object")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
seReference =readRDS(args$control)
seDisease=readRDS(args$treanment)
se=cbind(seReference,seDisease)

lsiReference=calcLSI(assay(seReference), nComponents = 50, binarize = TRUE, nFeatures = NULL)
lsiDisease=calcLSI(assay(seDisease), nComponents = 50, binarize = TRUE, nFeatures = NULL)

#Identify Promoter Overlapping Peaks +/- TSS 500 Bp
tssmm10 <- TxDb.Mmusculus.UCSC.mm10.knownGene %>% transcripts %>% resize(1,"start") %>% unique
strand(tssmm10) <- "*"
promoterPeaks <- subjectHits(findOverlaps(resize(tssmm10, 500 * 2 + 1), rowRanges(se), ignore.strand=TRUE))

#Input Parameters
input_knn <- 25
scaleTo <- 10000
nMax <- 500

#LSI-SVD
svdReference <- as.data.frame(lsiReference$matSVD)
#svdDisease <- as.data.frame(as.matrix(lsiProjection$matSVD))
svdDisease<-as.data.frame(lsiDisease$matSVD)
#Differential Seed
set.seed(1)

#Cells that we are testing of disease
plotDF=readRDS("CallPeaks_Matrix/plotDF.rds")
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
matHealthy <- assay(se[,idxReference])
matDisease <- assay(se[,idxDisease])

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

plotDir <- args$outdir
pdf(file.path(plotDir,"Differential-MA-Plot.pdf"), width = 8, height = 6, useDingbats = FALSE)

ggplot(plotDiff, aes(log2Mean,log2FC,color=type)) +
    ggrastr::geom_point_rast(size=0.5) +
    theme_bw() +
    xlab("log2 Mean") +
    ylab("log2 Fold Change") +
    scale_color_manual(values=c("not-differential"="lightgrey", "do-regulated"="dodgerblue3", "up-regulated"="firebrick3"))

dev.off()

#Save Output
readr::write_tsv(dfTT, file.path(plotDir,"Differential-Results.tsv"))
