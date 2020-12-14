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

message("INFO : Loading dataset")
se=readRDS(args$se)
counts=assay(se)
metadata=as.data.frame(colData(se))

message("INFO : Create Seurat Object")
seurat=CreateSeuratObject(counts,project="scATAC-GeneScore",meta.data=metadata)
message("INFO : Normalize")
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)

message("INFO : Save ")
saveRDS(seurat,file.path(args$outdir,"seurat.rds"))

message("INFO : Done!")

