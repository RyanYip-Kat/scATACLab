##这一步目的是先将RNA和ATAC变成RSE
#该脚本必须要ArchR1.0.0下进行
#library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(optparse)
library(yaml)
library(Rcpp)
library(Seurat)
library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(viridis)
set.seed(1)

####################################################
#Functions
####################################################
source('/home/shi/Work/ATAC/imd/Runcode-temp/integra-output/p2g_source.R')


####################################################
#Input Data
####################################################





parser$add_argument("--archr",
                    type="character",
                    default=NULL,
                    help="the path of project saved")

parser$add_argument("--seuratSE",
                    type="character",
                    default=NULL,
                    help="the path of project saved")

parser$add_argument("--outdir",
                    type="character",
                    default="")

parser$add_argument("--name",
                    type="character",
                    default="pbmc",
                    help="the dataset  to be used")       

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}


#rna需要先用 /home/shi/Work/ATAC/imd/Seurat_tohome/shi/Work/ATAC/imd/脚本换成SE格式
addArchRThreads(threads = 1)
rna=readRDS(args$seuratSE)
atac=loadArchRProject(args$archr)
name=args$name
seB=rna
proj=atac

#Window around TSS to be called promoter
tssWindow <- 2500

#Flanking distance from TSS in KB for Co-Accessibility
gtf_file='/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf'

tssWindow <- 2500
genes <- getGeneGTF(gtf_file) %>% resize(1,"start") %>% resize(tssWindow * 2 + 1, "center")
##seRNA TO RseRNA
names(genes) <- genes$gene_name
features=rownames(seB)
features=features[features%in%mcols(genes)$gene_name]
seB=seB[features,]
rowranges = genes[rownames(seB),]
rowRanges(seB)=rowranges
saveRDS(seB,file.path(args$outdir,paste0(args$name,"_RangedSE_RNA.rds")))

#这里是转atac
peaks=getPeakSet(proj)
seA=getMatrixFromProject(proj,'PeakMatrix')
rowRanges(seA)=peaks
rgs=as.data.frame(ranges(peaks))
sqs=as.data.frame(seqnames(peaks))
peakset=cbind(sqs,rgs)
peak_set=paste(peakset$value,peakset$start,peakset$end,sep="_")
rownames(seA)=peak_set
saveRDS(seA,file.path(args$outdir,paste0(args$name,"_RangedSE_ATAC.rds")))


