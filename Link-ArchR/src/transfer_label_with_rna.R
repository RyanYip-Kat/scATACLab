library(ggplot2)
library(RColorBrewer)
library(argparse)
library(Seurat)
library(Signac)
library(Cairo)
library(viridis)
library(future)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--rna",
                    type="character",
                    default=NULL)


parser$add_argument("--atac",
                    type="character",
                    default=NULL)


parser$add_argument("--ref_label",
                    type="character",
                    default=NULL)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

plan("multiprocess", workers = 12)

print("### Loading data")
rna=readRDS(args$rna)
atac=readRDS(args$atac)
DefaultAssay(atac) <- 'RNA'

rna <- FindVariableFeatures(
  object = rna,
  nfeatures = 5000
)

print("### Transfer anchors")
transfer.anchors <- FindTransferAnchors(
  reference =rna,
  query = atac,
  reduction = 'cca',
  dims = 1:40
)

print("### Transfer label")
metadata=rna$meta.data
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = metadata[[args$ref_label]],
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

atac <- AddMetaData(object = atac, metadata = predicted.labels)

print("### Save")
saveRDS(atac,file.path(args$outdir,"seurat.rds"))

