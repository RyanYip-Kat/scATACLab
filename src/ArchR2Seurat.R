library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
library(edgeR)
library(matrixStats)
library(SummarizedExperiment)
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(harmony)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

parser$add_argument("--nTop",
                    type="integer",
                    default=50000,
                    help="number of top var peaks")


parser$add_argument("--binarize",
		    action="store_true")


parser$add_argument("--peaktype",
		    nargs="+",
                    type="character",
                    default=NULL,
		    help="the peakType in PeakMatrix mcols for subset")

parser$add_argument("--batch_key",
                    type="character",
                    default=NULL,
                    help="barcode file use for subsetting")

parser$add_argument("--batch_correct",
                    action='store_true', default=FALSE)


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

getReducedDimsNames<-function(project){
        return(names(project@reducedDims))
}


getEmbeddingsNames<-function(project){
        return(names(project@embeddings))
}


print(paste0("### Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("### Loading Project")
projHeme<-loadArchRProject(args$project)


################# Peak Matrix
avaliable_matrixs=getAvailableMatrices(projHeme)
stopifnot("PeakMatrix"%in%avaliable_matrixs)
se=getMatrixFromProject(projHeme,useMatrix="PeakMatrix")

peaks=getPeakSet(projHeme)
rowRanges(se)=peaks
# mcols(peaks)
rowData(se)=cbind(rowData(se),mcols(peaks))
rgs=as.data.frame(ranges(peaks))
sqs=as.data.frame(seqnames(peaks))
peakset=cbind(sqs,rgs)
peak_set=paste(peakset$value,peakset$start,peakset$end,sep="-")
rownames(se)=peak_set
print(head(rowRanges(se)))

if(!is.null(args$peaktype)){
	se=subset(se,peakType%in%args$peaktype)
}

print("### Add ReducedDims")
ReducedDims_Names=c("IterativeLSI","Harmony")
for(name in ReducedDims_Names){
	slot=getReducedDims(projHeme,reducedDims=name,returnMatrix=F)
	metadata(se)[[name]]=slot
}

print(names(metadata(se)))
print("### Add Embedding")
mat=as.matrix(getEmbedding(ArchRProj =projHeme, embedding ="UMAP", returnDF = TRUE))
#Add UMAP coordinates to column data in summarized experiment
colData(se)$UMAP1 <-mat[,1]
colData(se)$UMAP2 <-mat[,2]

mat=as.matrix(getEmbedding(ArchRProj =projHeme, embedding ="TSNE", returnDF = TRUE))
colData(se)$TSNE1 <-mat[,1]
colData(se)$TSNE2 <-mat[,2]


###################
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}
###############
#Make Pseudo Bulk Library
mat=assay(se)
if(args$binarize){
	message("Making PseudoBulk...")
	mat@x[mat@x > 0] <- 1 #binarize
}

cluster=colData(se)$Clusters
clusterSums <- groupSums(mat = mat, groups = cluster, sparse = TRUE) #Group Sums
logMat <- edgeR::cpm(clusterSums, log = TRUE, prior.count = 3) #log CPM matrix
varPeaks <- head(order(matrixStats::rowVars(logMat), decreasing = TRUE), args$nTop) #Top variable peaks
metadata(se)$variablePeaks <- varPeaks


print("### Save Matrix")
saveRDS(se,file.path(args$outdir,"PeakMatrix_Summarized-Experiment.rds"))

message("### Create ATAC Seurat Object")
meta=as.data.frame(colData(se))

print("### Create Seurat")
seurat=CreateSeuratObject(mat,assay="peaks")
print(dim(seurat))
seurat=AddMetaData(seurat,metadata=meta)
print("### Add Eembeddings ")
#ReducedDims_Names=c("IterativeLSI","Harmony")
ReducedDims_Names=getReducedDimsNames(projHeme)
for(name in ReducedDims_Names){
        print(paste0("#### Add ",name," into seurat reduction"))
        mat=getReducedDims(projHeme,reducedDims=name,returnMatrix=T)
        slot=paste("ar",str_to_lower(name),sep="")
        seurat[[slot]]<-CreateDimReducObject(embeddings =mat,
                                  key = paste0(slot,"_"),
                                  assay = DefaultAssay(seurat))
}

#Embedding_Names=c("UMAP","TSNE")
Embedding_Names=getEmbeddingsNames(projHeme)
for(name in Embedding_Names){
        print(paste0("#### Add ",name," into seurat reduction"))
        mat=as.matrix(getEmbedding(ArchRProj =projHeme, embedding =name, returnDF = TRUE))
        slot=paste("ar",str_to_lower(name),sep="")
        colnames(mat)=paste(paste0(slot,"_"),1:ncol(mat),sep = "")
        seurat[[slot]]<-CreateDimReducObject(embeddings =mat,
                                  key = paste0(slot,"_"),
                                  assay = DefaultAssay(seurat))
}

print("### Run Eembedding and Clusters with Seurat")
seurat <- RunTFIDF(seurat)
seurat <- FindTopFeatures(seurat, min.cutoff = 'q0')
seurat <- RunSVD(seurat)

use_rep="lsi"
nDims=2:30
if(args$batch_correct){
        vars=ifelse(!is.null(args$batch_key),args$batch_key,"Sample")
        stopifnot(vars%in%colnames(seurat@meta.data))
        print(paste0("### Run Harmony with :",vars))
        #seurat=RunHarmony(seurat, group.by.vars =vars,
        # 			  reduction = "lsi",dims.use=2:30,assay.use=DefaultAssay(seurat))
	metadata=seurat@meta.data
        m=Embeddings(seurat,"lsi")[,nDims]
	harmony_embeddings=HarmonyMatrix(m, metadata, vars) #v3
	rownames(harmony_embeddings)=colnames(seurat)
        colnames(harmony_embeddings)=paste("Harmony",1:ncol(harmony_embeddings),sep="_")
        mat=as.matrix(harmony_embeddings)
        colnames(mat)<-paste("Harmony_",1:ncol(mat),sep = "")
        seurat[["harmony"]]<-CreateDimReducObject(embeddings =mat,
                                                  key = "Harmony_",
                                                  assay = DefaultAssay(seurat))

	nDims=1:20
	use_rep="harmony"
}

seurat <- RunUMAP(object = seurat, reduction =use_rep, dims = nDims)
seurat <- FindNeighbors(object = seurat, reduction = use_rep, dims = nDims)
seurat <- FindClusters(object = seurat, verbose = FALSE, algorithm = 3)

#######################    Gene Score
print("### Get Score GeneMatrix")
GS=getMatrixFromProject(projHeme,useMatrix = "GeneScoreMatrix")
features=getFeatures(projHeme,useMatrix = "GeneScoreMatrix")
rownames(GS)=features

score=assay(GS)

print("### Create Seurat")
seurat[['RNA']] <- CreateAssayObject(counts = score)
seurat <- NormalizeData(
  object = seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_RNA)
)
DefaultAssay(seurat)="RNA"
seurat <- FindVariableFeatures(seurat, selection.method = "vst",
                            nfeatures = 2000,verbose = FALSE)


print("### Save Seurat")
saveRDS(seurat,file.path(args$outdir,"ArchR_Seurat.rds"))

