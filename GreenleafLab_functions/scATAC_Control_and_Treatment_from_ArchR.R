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

parser$add_argument("--split_by",
                    type="character",
                    default="Status",
		    help="which column as the control and treatment")

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
se=readRDS(args$combinedSE)
cellColData=as.data.frame(colData(se))
split_by=cellColData[[args$split_by]]
#Plot Info

cells=ifelse(split_by==args$control,"reference","disease")
print(table(cells))
clusters=cellColData$Clusters
splitCells <- split(cells,clusters)
df <- data.frame(
    clusters = names(splitCells),
    proportion = unlist(lapply(seq_along(splitCells), function(x) sum(splitCells[[x]]=="disease") / length(splitCells[[x]])))
    )


print("### Plot UMAP Data Frame")
plotDF <- data.frame(X1=cellColData$UMAP1,X2=cellColData$UMAP2)
rownames(plotDF) <- rownames(cellColData)
plotDF$type <- cells
plotDF$clusters <- clusters
plotDF$classification <- 0
#If disease cells are clustered with healthy cluster (proportion > 0.9) we will classify these as healthy-like
plotDF$classification[plotDF$type == "disease" & plotDF$clusters %in% paste0(df$clusters[df[,2] > args$ratio])] <- 1
plotDF$classification[plotDF$type == "disease"] <- plotDF$classification[plotDF$type == "disease"] + 1
plotDF <- plotDF[order(plotDF$classification), ]

#Formal Classification
plotDF$classificationSTR <- "reference"
plotDF$classificationSTR[plotDF$classification==1] <- "healthy-like"
plotDF$classificationSTR[plotDF$classification==2] <- "disease-like"

#Plot PDFs
plotDir <- args$outdir
pdf(file.path(plotDir,"Classification-UMAP.pdf"), width = 12, height = 12, useDingbats = FALSE)

ggplot(plotDF, aes(X1,X2,color=classificationSTR)) +
    geom_point() +
    theme_bw() +
    xlab("UMAP Dimension 1") +
    ylab("UMAP Dimension 2") +
    scale_color_manual(values=c("reference"="lightgrey","healthy-like"="dodgerblue3","disease-like"="firebrick3"))

dev.off()

saveRDS(plotDF,file.path(args$outdir,"plotDF.rds"))

# Ge promoterPeakst
tss <- TxDb.Mmusculus.UCSC.mm10.knownGene %>% transcripts %>% resize(1,"start") %>% unique
strand(tss) <- "*"
promoterPeaks <- subjectHits(findOverlaps(resize(tss, 500 * 2 + 1), rowRanges(se), ignore.strand=TRUE))

saveRDS(promoterPeaks,file.path(args$outdir,"promoterPeaks.rds"))
