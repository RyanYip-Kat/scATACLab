library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(argparse)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--path",
                    type="character",
                    default=NULL,
		    help="The path from cellranger-atac count or aggr")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")


parser$add_argument("--archr_barcode",
                    type="character",
                    default=NULL,
		    help="the barcode from ArchR handled barcode")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

fragment=file.path(args$path,"outs/fragments.tsv.gz")
singlecell=file.path(args$path,"outs/singlecell.csv")
matrix_h5=file.path(args$path,"outs/filtered_peak_bc_matrix.h5")
print("### Pre-processing workflow and loading dataset")
counts <- Read10X_h5(matrix_h5)
metadata <- read.csv(
  file = singlecell,
  header = TRUE,
  row.names = 1
)

seurat_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome =args$genome,
  fragments = fragment,
  min.cells = 1
)

seurat <- CreateSeuratObject(
  counts = seurat_assay,
  assay = 'peaks',
  project = 'ATAC',
  min.cells=200,
  meta.data = metadata
)

print(paste0("### Annotation with ",args$genome))
# extract gene annotations from EnsDb
if(args$genome=="hg38"){
	annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
}else{
	annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
}

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- args$genome
# add the gene information to the object
Annotation(seurat) <- annotations
print(paste0("### The original object shape : [",ncol(seurat),",",nrow(seurat),"]"))

print("### Subset with ArchR QC barcodes")
preData=read.csv(args$archr_barcode,sep=",")
preMeta.col.name=colnames(preData)[2]
colnames(preData)=c("barcode","label")
rownames(preData)=preData$barcode
preMeta=subset(preData,select=label)

seurat=subset(seurat,cells=rownames(preMeta))
seurat<-AddMetaData(seurat,metadata=preMeta,col.name=preMeta.col.name)
print(paste0("### After subset from ArchR QC's object shape : [",ncol(seurat),",",nrow(seurat),"]"))


print("### Saving")
saveRDS(seurat,file.path(args$outdir,"signac.rds"))

