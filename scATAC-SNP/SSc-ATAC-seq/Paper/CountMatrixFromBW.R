library(derfinder)
library(ArchR)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(GenomeInfoDb)
library(rtracklayer)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(argparse)
library(yaml)
library(Rcpp)

source("/home/ye/Work/BioAligment/SNP/Shi/SSc-ATAC-seq/Paper/functions.R")
makedir=function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

########################################   callpeak paramenters
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--meta",
                    type="character",
                    default=NULL,
                    help="*.bw file in meta csv file")


args <- parser$parse_args()

outDir=args$outdir
makedir(outDir)

##################
DF=read.csv(args$meta,header=FALSE,sep=",",stringsAsFactors=FALSE)
bwFiles=as.character(DF$V2)
Names=as.character(DF$V1)
names(bwFiles)=Names

# change to UCSC style since the data was mapped to hg38
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
annotations=filterChrGR(annotations,remove=c("chrM"))
exons=annotations[mcols(annotations)$type=="exon"]

#txdb <- keepSeqlevels(TxDb.Hsapiens.UCSC.hg38.knownGene, paste("chr",c(1:22,"X","Y"),sep=""))
#exons <- exonsBy(txdb, by = 'gene')

message("INFO : Compute Matrix Counts")
bw <- BigWigFileList(bwFiles)
counts <- matrix(NA, nrow = length(exons), ncol = length(bwFiles))
colnames(counts) <- names(bw)
for(i in seq_len(length(bw))) {
    cat(sprintf("%d of %d",i,length(bw)))
    coverage <- import(bw[[i]], as = 'RleList')$chr21
    counts[, i] <- sum(Views(coverage, ranges(exons)))
}


message("INFO : Save ...")
saveRDS(counts,file.path(outDir,"bwCounts.rds"))

message("INFO : Done!")
