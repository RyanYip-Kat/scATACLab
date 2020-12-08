library(Seurat)
library(Signac)
library(argparse)
library(stringr)
library(future)

library(SeuratWrappers)
library(monocle3)
library(Matrix)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default="",
                    help="the path of seurat has been created")

parser$add_argument("--batch",
                    type="character",
                    default="Status",
                    help="which column to use subset")

#parser$add_argument("--invert",
#                    action='store_true', default=FALSE)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


print("### Loading Dataset")
seurat<-readRDS(args$seurat)
seurat <- FindTopFeatures(seurat)
seurat<-subset(seurat,features=VariableFeatures(seurat))

stopifnot(args$batch%in%colnames(seurat@meta.data))
counts<-GetAssayData(seurat,"counts")
print("#### gene meta data")
pd<-seurat@meta.data
rownames(pd)<-colnames(counts)
fd <- data.frame(gene_short_name = row.names(counts), row.names = row.names(counts))
print("#### new cell data set")
cds<-new_cell_data_set(counts,cell_metadata=pd,gene_metadata=fd)
cds<-detect_genes(cds)

print("### preprocess ")
cds <- preprocess_cds(cds,
                      num_dim = 50,
                      method="LSI",
                      norm_method="log")

print("### Align")
cds <- align_cds(cds,
                 preprocess_method="LSI",
                 alignment_k=20,
                 residual_model_formula_str="~nCount_peaks+BlacklistRatio+TSSEnrichment",
                 alignment_group=args$batch)


print("### reduce dimension")
cds <- reduce_dimension(cds,reduction_method="tSNE",preprocess_method="Aligned",cores=8)
cds <- reduce_dimension(cds,reduction_method="UMAP",preprocess_method="Aligned",cores=8)

print("### cluster")
cds<-cluster_cells(cds,
                   reduction_method="UMAP",
                   k=20,
                   cluster_method="leiden",
                   partition_qval=0.05)

cds<-cluster_cells(cds,
                   reduction_method="tSNE",
                   k=20,
                   cluster_method="leiden",
                   partition_qval=0.05)

print("### learn graph")
cds<-learn_graph(cds,
                 use_partition=TRUE,
                 close_loop=TRUE)

print("### Save monocle")
saveRDS(cds,file.path(args$outdir,"monocle.rds"))

#seurat=SCTransform(seurat,vars.to.regress="nFeature_RNA",verbose = FALSE)
#########################
print("### Create ReducedDim from monocle and add clusters")
#tSNE_clusters<-clusters(cds,reduction_method="tSNE")
UMAP_clusters<-clusters(cds,reduction_method="UMAP")

#tSNE_partitions<-partitions(cds,reduction_method="tSNE")
UMAP_partitions<-partitions(cds,reduction_method="UMAP")

monocle_meta<-data.frame(
                         "UMAP_clusters"=UMAP_clusters,
                         "UMAP_partitions"=UMAP_partitions,
                         row.names=names(UMAP_clusters))

print("### Add MetaData")
seurat<-AddMetaData(seurat,metadata=monocle_meta)

print("### Add reducedDims")
print("#### Add UMAP")
mat<-reducedDims(cds)[["UMAP"]]
colnames(mat)<-paste("UMAP_",1:ncol(mat),sep = "")
seurat[["umap"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "UMAP_",
                                  assay = DefaultAssay(seurat))

print("#### Add pca")
mat<-reducedDims(cds)[["LSI"]]
colnames(mat)<-paste("LSI_",1:ncol(mat),sep = "")
seurat[["lsi"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "LSI_",
                                  assay = DefaultAssay(seurat))

print("#### Add Aligned")
mat<-reducedDims(cds)[["Aligned"]]
colnames(mat)<-paste("Aligned_",1:ncol(mat),sep = "")
seurat[["aligned"]]<-CreateDimReducObject(embeddings =mat,
                                  key = "Aligned_",
                                  assay = DefaultAssay(seurat))

saveRDS(seurat,file.path(args$outdir,"seurat_monocle3.rds"))


