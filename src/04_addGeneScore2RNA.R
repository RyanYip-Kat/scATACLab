library(Seurat)
library(Signac)
library(argparse)
library(stringr)
library(future)

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

plan("multiprocess", workers = 16)
outDir=args$outdir
Genome=args$genome

makedir(outDir)
if(Genome=="hg38"){
        library(EnsDb.Hsapiens.v86)
        annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
}else{
        library(EnsDb.Mmusculus.v79)
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
}

# change to UCSC style since the data was mapped to hg19
cat(sprintf("INFO : Annotation with [ %s ]\n",Genome))
seqlevelsStyle(annotations) <- 'UCSC'

############################## Loading dataset
message("INFO : Loading dataset")
object=readRDS(args$seurat)
Annotation(object) <- annotations
##############################  Create a gene activity matrix 
message("INFO : Get GeneActivity")
gene.activities <- GeneActivity(object)
# add the gene activity matrix to the Seurat object as a new assay and normalize it

message("INFO : Add RNA Slot with  GeneActivity ...")
object[['RNA']] <- CreateAssayObject(counts = gene.activities)
object <- NormalizeData(
  object = object,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(object$nCount_RNA)
)

message("INFO : Save Object ...")
saveRDS(object,args$seurat)
