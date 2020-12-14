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
parser <- ArgumentParser(description='Program to convert Gene Activity Object into Seurat Object')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--se",
                    type="character",
                    default=NULL,
                    help="Cicero Gene Activity object")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
outDir=args$outdir
makedir(outDir)
annoFile="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"

message("INFO : Loading dataset")
se=readRDS(args$se)
counts=assay(se)
rownames(counts)=str_replace_all(rownames(counts),"_","-")
seurat=CreateSeuratObject(counts,min.cells=1)
counts=GetAssayData(seurat,"counts")
metadata=as.data.frame(colData(se))

activity.matrix <- CreateGeneActivityMatrix(peak.matrix = counts, annotation.file = annoFile,
    seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE,keep.sparse=FALSE)

message("INFO : Add RNA Slot with  Gene Activity ...")
seurat <- CreateSeuratObject(counts = log2(10^6*activity.matrix+1),assay="RNA",meta.data=metadata)
seurat <- NormalizeData(
  object = seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_RNA)
)


message("INFO : Save ")
saveRDS(seurat,file.path(args$outdir,"seurat.rds"))

message("INFO : Done!")

