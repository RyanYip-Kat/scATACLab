library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(argparse)
library(harmony)
library(sctransform)
library(SeuratWrappers)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--object",
                    type="character",
                    default=NULL,
		    help="The seurat object get from preprocess with ArchR QC Barcode")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--batch_key",
                    type="character",
                    default="sample",
		    help="the batch key for batch correct")

parser$add_argument("--batch_correct",
		    action="store_true",
		    default=FALSE)

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print(paste0("### Loading  object from :",args$object))
seurat=readRDS(args$object)

print("### Normalization and linear dimensional reduction")
#seurat <- SCTransform(seurat, assay="peaks",verbose = TRUE,variable.features.n=10000)
seurat <- RunTFIDF(seurat)
seurat <- FindTopFeatures(seurat, min.cutoff = 'q0')
seurat <- RunSVD(object = seurat)

use_rep="lsi"
if(args$batch_correct){
	print(paste0("### Correct Batch effect with : ",args$batch_key," by Harmony"))
	metadata=seurat@meta.data
	stopifnot(args$batch_key%in%colnames(metadata))
	seurat <- SCTransform(seurat, assay="peaks",verbose = TRUE,variable.features.n=10000)
	seurat<- RunHarmony(seurat, group.by.vars =args$batch_key,reduction="lsi",assay="SCT",dims.use = 2:30)
	use_rep="harmony"

}

print("### Run UMAP") 
if(use_rep=="lsi"){
	dims.use=2:30
}else{
	dims.use=1:20
}
seurat <- RunUMAP(
  object = seurat,
  reduction =use_rep,
  dims = dims.use
)
print("### Find Neigbors")
seurat <- FindNeighbors(
  object = seurat,
  reduction = use_rep,
  dims = dims.use
)
print("### Find Clusters")
seurat <- FindClusters(
  object = seurat,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

print("### Create a gene activity matrix")
gene.activities <- GeneActivity(seurat)

# add the gene activity matrix to the Seurat object as a new assay
seurat[['RNA']] <- CreateAssayObject(counts = gene.activities)
seurat <- NormalizeData(
  object = seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_RNA)
)

print("### Save")
saveRDS(seurat,file.path(args$outdir,"signac_embedding_cluster.rds"))
