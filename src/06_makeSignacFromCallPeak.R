library(stringr)
library(DropletUtils)
library(ggplot2)
library(Seurat)
library(Signac)
library(argparse)
library(GenomicRanges)
library(future)
#############################
source("/home/ye/Work/R/scATAC/Signac/src/helper_functions.R")
parser <- ArgumentParser(description='Program to make Signac Object')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="seurat Object")

parser$add_argument("--peak",
                    type="character",
                    default=NULL,
                    help="result from callpeaks")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")

args <- parser$parse_args()

###############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}


###############################  Configure
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
genome(annotations) <- Genome


################################  loading dataset
seurat=readRDS(args$seurat)
peaks=readRDS(args$peak)

peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

message("INFO : get callpeak counts ...")
macs2_counts=FeatureMatrix(fragments=Fragments(seurat),features=peaks,cells = colnames(seurat))

message("INFO : Create Object with callpeak counts ...")
metadata=seurat@meta.data
peaks.assay=CreateChromatinAssay(counts=macs2_counts,
                                  fragments=Fragments(seurat),
                                  genome =Genome,
                                  sep = c("-", "-"),
                                  annotation=annotations,
                                  validate.fragments = FALSE)
seurat.peaks <- CreateSeuratObject(peaks.assay, assay = "peaks",meta.data = metadata)
Annotation(seurat.peaks) <- annotations
seurat.peaks[["ATAC"]]=seurat[["ATAC"]]
message("INFO : Save Object ...")
saveRDS(seurat.peaks,file.path(outDir,"seurat-peaks.rds"))
message("INFO : Done!")
