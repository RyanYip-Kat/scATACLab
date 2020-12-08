library(Seurat)
library(Signac)
library(argparse)
library(stringr)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")


parser$add_argument("--outdir",
                    type="character",
                    default="./Results")


parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")


args <- parser$parse_args()

############################### Configure
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

outDir=args$outdir
Genome=args$genome

makedir(outDir)
if(Genome=="hg38"){
        #library(EnsDb.Hsapiens.v86)
        #annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
        annoFile="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
}else{
        #library(EnsDb.Mmusculus.v79)
        #annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
        annoFile="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-mm10-3.0.0/genes/genes.gtf"
}

# change to UCSC style since the data was mapped to hg19
#cat(sprintf("INFO : Annotation with [ %s ]\n",Genome))
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- Genome

############################## Loading dataset
message("INFO : Loading dataset")
seurat=readRDS(args$seurat)
counts=GetAssayData(seurat,"counts")
message("INFO : Get GeneActivity")
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = counts, annotation.file = annoFile,
    seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE,keep.sparse=FALSE)

message("INFO : Add RNA Slot with  Gene Activity ...")
seurat[['RNA']] <- CreateAssayObject(counts = activity.matrix)
seurat <- NormalizeData(
  object = seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat$nCount_RNA)
)


message("INFO : Save Object ...")
saveRDS(seurat,args$seurat)
message("INFO : Done!")
