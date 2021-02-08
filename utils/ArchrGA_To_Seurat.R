library(argparse)
library(stringr)
library(Seurat)
library(Signac)
library(ArchR)
library(Matrix)

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

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("### Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("### Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)

#matrixs=c("PeakMatrix","GeneScoreMatrix")
#avaliable_matrixs=getAvailableMatrices(projHeme)
#print(avaliable_matrixs)

print("### Get GeneMatrix")
RNA=getMatrixFromProject(projHeme,useMatrix = "GeneScoreMatrix")
score_file=file.path(getOutputDirectory(projHeme),"ArchR-GeneMatrix-Project.rds")
if(!file.exists(score_file)){
	print("### Save GeneScore Matrix")
	saveRDS(RNA,score_file)
}


features=getFeatures(projHeme,useMatrix = "GeneScoreMatrix")
rownames(RNA)=features

meta=as.data.frame(colData(RNA))
score=assay(RNA)

print("### Create Seurat")
seurat=CreateSeuratObject(score)
seurat=AddMetaData(seurat,metadata=meta)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst",
                            nfeatures = 5000,verbose = FALSE)

print("### Add Eembeddings ")
ReducedDims_Names=c("IterativeLSI","Harmony")
for(name in ReducedDims_Names){
	print(paste0("#### Add ",name," into seurat reduction"))
	mat=getReducedDims(projHeme,reducedDims=name,returnMatrix=T)
	slot=str_to_lower(name)
	seurat[[slot]]<-CreateDimReducObject(embeddings =mat,
                                  key = paste0(slot,"_"),
                                  assay = DefaultAssay(seurat))
}

Embedding_Names=c("UMAP","TSNE")
for(name in Embedding_Names){
	print(paste0("#### Add ",name," into seurat reduction"))
	mat=as.matrix(getEmbedding(ArchRProj =projHeme, embedding =name, returnDF = TRUE))
	slot=str_to_lower(name)
	colnames(mat)=paste(paste0(slot,"_"),1:ncol(mat),sep = "")
	seurat[[slot]]<-CreateDimReducObject(embeddings =mat,
                                  key = paste0(slot,"_"),
                                  assay = DefaultAssay(seurat))
}

print("### Save Seurat")
saveRDS(seurat,file.path(args$outdir,"ArchR_GA_Seurat.rds"))
