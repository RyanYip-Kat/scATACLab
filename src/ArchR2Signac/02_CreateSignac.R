library(stringr)
library(DropletUtils)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(argparse)
library(GenomicRanges)
library(future)


source("/home/ye/Work/BioAligment/SNP/Shi/scripts/help_functions.R")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--path",
                    type="character",
                    default=NULL,
                    help="path from  ArchR prepare for Signac")


parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")

args <- parser$parse_args()

################################
bgzip="/home/ye/anaconda3/envs/BulkBio/bin/bgzip"
tabix="/home/ye/anaconda3/envs/BulkBio/bin/tabix"

outDir=args$path
path=args$path
Genome=args$genome

################################# configure

fileNames=list.files(path)
cat(sprintf("INFO : There are [ %d ] files in [ %s ] \n",length(fileNames),path))
message("INFO : loading datasets for Creating Signac Object ...")
assertthat::assert_that("seurat_list.rds"%in%fileNames)
assertthat::assert_that("fragment_assay_list.rds"%in%fileNames)
assertthat::assert_that("meta_list.rds"%in%fileNames)
assertthat::assert_that("gA.rds"%in%fileNames)
assertthat::assert_that("fragment_list.rds"%in%fileNames)
assertthat::assert_that("cell_list.rds"%in%fileNames)
assertthat::assert_that("ReducedDims"%in%fileNames)
assertthat::assert_that("Embeddings"%in%fileNames)

if(Genome=="hg38"){
        library(EnsDb.Hsapiens.v86)
        annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
        annoFile="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
}else{
        library(EnsDb.Mmusculus.v79)
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
        annoFile="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-mm10-3.0.0/genes/genes.gtf"
}
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- Genome

################################## load datasets
message("INFO : loading datasets for Creating Signac Object ...")
cat(sprintf("INFO : Loading Seurat list from --- [ %s ] ... \n",file.path(path,"seurat_list.rds")))
seurat_list=readRDS(file.path(path,"seurat_list.rds"))

cat(sprintf("INFO : Loading fragment assay list from --- [ %s ] ... \n",file.path(path,"fragment_assay_list.rds")))
fragment_assay_list=readRDS(file.path(path,"fragment_assay_list.rds"))

cat(sprintf("INFO : Loading ArchR gA from --- [ %s ] ... \n",file.path(path,"gA.rds")))
score=readRDS(file.path(path,"gA.rds"))

cat(sprintf("INFO : Loading metadata list from --- [ %s ] ... \n",file.path(path,"meta_list.rds")))
meta_list=readRDS(file.path(path,"meta_list.rds"))

###################################
message("INFO : Merge Count lists ...")
#seurat=Reduce(merge,seurat_list)
count_list=lapply(seurat_list,function(seurat){return(GetAssayData(seurat,"counts",assay="ATAC"))})
counts=do.call(cbind,count_list)
##################################### Create Combined Object
message("INFO : Make Combined Object ...")
#counts=GetAssayData(seurat,"counts",assay="ATAC")
fragments=fragment_assay_list
metadata=do.call(rbind,meta_list)
#metadata=seurat@meta.data

message("INFO : CreateChromatinAssay ...")
peaks.assay=CreateChromatinAssay(counts=counts,
                                  fragments=fragments,
                                  genome =args$genome,
                                  sep = c("-", "-"),
                                  annotation=annotations,
                                  validate.fragments = FALSE)
seurat <- CreateSeuratObject(peaks.assay, assay = "ATAC",meta.data = metadata)

message("INFO : gA RNA Assay Object ...")
archrGA <- CreateAssayObject(counts =score)
seurat[['archrGA']] <- archrGA
seurat <- NormalizeData(object = seurat,
                        assay = 'archrGA',
                        normalization.method = 'LogNormalize')


############################### Add ArchR reduced slots
ReducedDimDir=file.path(path,"ReducedDims")
EmbeddingDir=file.path(path,"Embeddings")

embeds=list.files(EmbeddingDir)
reuds=list.files(ReducedDimDir)
prefix="ar"

message("INFO : Add ArchR ReduceDims into Signac Object")
for(name in embeds){
	slotName=str_split(name,"\\.")[[1]][1]
	cat(sprintf("INFO : add Slot [ %s ]\n",slotName))
        mat=readRDS(file.path(EmbeddingDir,name))

        slot=paste0(prefix,str_to_lower(slotName))
        #colnames(mat)=paste(paste0(slot,"_"),1:ncol(mat),sep = "")
        seurat[[slot]]<-CreateDimReducObject(embeddings =mat,
                                  key = paste0(slot,"_"),
                                  assay = "ATAC")
}

for(name in reuds ){
	slotName=str_split(name,"\\.")[[1]][1]
	cat(sprintf("INFO : add Slot [ %s ]\n",slotName))
	mat=readRDS(file.path(ReducedDimDir,name))

        slot=paste0(prefix,str_to_lower(slotName))
        seurat[[slot]]<-CreateDimReducObject(embeddings =mat,
                                  key = paste0(slot,"_"),
                                  assay = "ATAC")
}
message("INFO : add Slot Done!")
message("INFO : Save Combined Object")
saveRDS(seurat,file.path(outDir,"combined_seurat.rds"))
message("INFO : Done!")
