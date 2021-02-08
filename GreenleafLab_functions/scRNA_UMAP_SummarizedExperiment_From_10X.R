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
library(Seurat)
set.seed(1)

functions="/home/ye/Work/R/scATAC/ArchR/GreenleafLab_functions/functions.R"
source(functions)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="./output")

parser$add_argument("--crg",
                    type="character",
                    default="cellranger count or aggr path")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
                    help="which column to use subset")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset")

parser$add_argument("--invert",
                    action='store_true', default=FALSE)

parser$add_argument("--soupX",
                    action='store_true', default=FALSE)

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
path=args$crg
toc=Read10X(file.path(path,"filtered_feature_bc_matrix"))
gp_cluster<-read.csv(file.path(path,"analysis/clustering/graphclust/clusters.csv"))
rownames(gp_cluster)<-gp_cluster$Barcode

km_cluster<-read.csv(file.path(path,"analysis/clustering/kmeans_10_clusters/clusters.csv"))
rownames(km_cluster)<-km_cluster$Barcode

gp_cluster<-subset(gp_cluster,select=Cluster)
km_cluster<-subset(km_cluster,select=Cluster)

print(paste0("Size of toc  [ ",nrow(toc),",",ncol(toc)," ]"))
genes<-rownames(toc)
keep_genes<-genes[!str_detect(genes,"^mt-|^Rpl|^Rps")]
keep_genes<-keep_genes[!str_detect(keep_genes,"\\.")]
toc=toc[keep_genes,]

if(args$soupX){
	library(SoupX)
	print("### Use SoupX method")
	tod=Read10X(file.path(path,"raw_feature_bc_matrix"))
	tod=tod[keep_genes,]
	print(paste0("Size of tod  [ ",nrow(tod),",",ncol(tod)," ]"))
	sc = SoupChannel(tod, toc)
	sc = setClusters(sc, setNames(gp_cluster$Cluster, rownames(gp_cluster)))
	sc = autoEstCont(sc)
        toc = adjustCounts(sc)
}


##################
seurat<-CreateSeuratObject(counts= toc,
                       assay = "RNA",
                       project ="scRNA",
                       names.delim="_",
                       min.cells=0,
                       min.features=0)

seurat<-AddMetaData(seurat,metadata=gp_cluster,col.name="orig.Cluster")
seurat<-AddMetaData(seurat,metadata=km_cluster,col.name="km.Cluster")
cells<-colnames(seurat)

ident=unlist(lapply(cells,function(cell){return(str_split(cell,"-")[[1]][2])}))
seurat$ident<-ident

metadata=seurat@meta.data
if(!is.null(args$column) && !is.null(args$subset)){
	print(paste0("### Subset with : ",args$column))
	Idents(seurat)<-metadata[[args$column]]
        seurat<-subset(seurat,idents=args$subset,invert=args$invert)
        print(paste0("after subset,number of cells : ",ncol(seurat)))
}

##################
seurat<-NormalizeData(seurat,normalization.method = "LogNormalize",verbose = FALSE)
logcount=GetAssayData(seurat,"data")
count=GetAssayData(seurat,"counts")

metadata=seurat@meta.data
print("### Convert into SummarizedExperiment")
se=SummarizedExperiment(assays=list(counts=count,logcounts=logcount),colData=metadata)
#se=SummarizedExperiment(assays=list(counts=count),colData=metadata)

nPCs <- 1:25 #Number of PCs for clustering
nTop <- c(3000, 3000, 3000) #Choose a higher number of variable features
resolution <- c(0.2,0.6,1.0) #Clustering resolutions for Seurat SNN

#Optimize LSI Features
print("### Optimize LSI Features")
lsiObj <- scRNA_optimizeLSI(assay(se),
  resolution = resolution,
  pcsUse = nPCs,
  varFeatures = nTop)


metadata(se)$optimizeLSI <- lsiObj
metadata(se)$matSVD <- lsiObj[[length(lsiObj)]][[1]][[1]] #Last one
metadata(se)$variableGenes <- rownames(se)[lsiObj[[length(lsiObj)]]$varFeatures] #Variable genes

####################################################
#For Creating UMAP Start Here
####################################################
matSVD <- metadata(se)$matSVD
clusters <- colData(se)$Clusters

#Set Seed and perform UMAP on LSI-SVD Matrix
print("### Run UMAP model")
set.seed(1)
uwotUmap <- uwot::umap(
    matSVD, 
    n_neighbors = 30, 
    min_dist = 0.45, 
    metric = "euclidean", 
    n_threads = 8, 
    verbose = TRUE, 
    ret_nn = TRUE,
    ret_model = TRUE
    )
save_uwot(uwotUmap,"UMAP-model.uwot")
system(paste0("mv UMAP-model.uwot ",args$outdir))
print("### Save umap model")
#Add UMAP coordinates to column data in summarized experiment
colData(se)$UMAP1 <- uwotUmap[[1]][,1]
colData(se)$UMAP2 <- uwotUmap[[1]][,2]
metadata(se)$uwot_model=file.path(args$outdir,"UMAP-model.uwot")
print("### Save SummarizedExperiment object")
saveRDS(se,file.path(args$outdir,"scRNA_SummarizedExperiment_Expriment.rds"))

