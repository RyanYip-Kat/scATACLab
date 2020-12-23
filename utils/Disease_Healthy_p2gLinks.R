##这一步是实现p2g
#该脚本必须要ArchR0.9.5下进行
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

parser$add_argument("--disease",
                    type="character",
                    default=NULL,
                    help="the name of disease status")

parser$add_argument("--healthy",
                    type="character",
                    default='HC',
                    help="the name of healthy")

parser$add_argument("--outdir",
                    type="character",
                    default="")


parser$add_argument("--seA",
                    type="character",
                    default=NULL,
                    help="RseATAC")

parser$add_argument("--seB",
                    type="character",
                    default=NULL,
                    help="RseRNA")                    

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
if(!dir.exists(args$outdir)){
  dir.create(args$outdir,recursive=TRUE)
}

#######################################
addArchRThreads(threads = 1)
proj=loadArchRProject(args$archr)
Total_proj=proj
disease=args$disease
healthy=args$healthy
seA=readRDS(args$seA)
seB=readRDS(args$seB)
###Disease
idxPass <- which(proj$Status %in% disease)
cellsPass <- proj$cellNames[idxPass]
Disease_proj <- proj[cellsPass, ]
###Healthy
idxPass <- which(proj$Status %in% healthy)
cellsPass <- proj$cellNames[idxPass]
Healthy_proj <- proj[cellsPass, ]


print(Healthy_proj)
print(Disease_proj)
print(Total_proj)

Healthy_proj <- addPeak2GeneLinks(
  ArchRProj = Healthy_proj,
  reducedDims = "Harmony"
)

Disease_proj <- addPeak2GeneLinks(
  ArchRProj = Disease_proj,
  reducedDims = "Harmony"
)

Total_proj <- addPeak2GeneLinks(
  ArchRProj = Total_proj,
  reducedDims = "Harmony"
)

saveRDS(Healthy_proj,file.path(args$outdir,paste0(healthy,"_P2gLinksProj.rds")))
saveRDS(Disease_proj,file.path(args$outdir,paste0(disease,"_P2gLinksProj.rds")))
saveRDS(Total_proj,file.path(args$outdir,paste0("Total_P2gLinksProj.rds")))

o_Disease=metadata(Disease_proj@peakSet)$Peak2GeneLinks
o_Healthy=metadata(Healthy_proj@peakSet)$Peak2GeneLinks
o_Total=metadata(Total_proj@peakSet)$Peak2GeneLinks


colnames(o_Total)=c('idxTotalA','idxTotalB','CorrelationTotal','FDRTotal','VarTotalA','VarTotalB')

colnames(o_Disease)=c('idxDiseaseA','idxDiseaseB','CorrelationDisease','FDRDisease','VarDiseaseA','VarDiseaseB')

colnames(o_Healthy)=c('idxHealthyA','idxHealthyB','CorrelationHealthy','FDRHealthy','VarHealthyA','VarHealthyB')

o=cbind(o_Total,o_Disease,o_Healthy)


fixA <- "center"
fixB <- "start"
associationWindow <- 2 * 250*10^3 + 1 #+-250 Kb
corCutOff <- 0.2 #Pearson Correlation Cutoff
fdrCutOff <- 0.1 #FDR Cutoff
distCutOff <- 2500 #Min Dist to TSS
gtf_file='/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf'
gtf <- getGeneGTF(gtf_file)
tssRNA <- resize(gtf, 1, "start")
strand(tssRNA) <- "*"
peakLinks <- rowRanges(seA)[o[,1]]
geneLinks <- rowRanges(seB) %>% resize(1, "start") %>% {.[o[,2]]}
mcolsLinks <- data.frame(geneLinks)[,c("seqnames","start","strand","gene_name","gene_id","exonLength")]
colnames(mcolsLinks) <- c("gene_chr","gene_start","gene_strand","gene_name","gene_id","exonLength")
mcolsLinks <- cbind(mcolsLinks, data.frame(o))
mcolsLinks$nearestGene <- tssRNA$gene_name[subjectHits(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))]
mcolsLinks$nearestGeneDistance <- mcols(distanceToNearest(peakLinks, tssRNA, ignore.strand=TRUE))$distance
mcols(peakLinks) <- mcolsLinks
peakLinks$peakName <- paste(seqnames(peakLinks), start(peakLinks), end(peakLinks), sep = "_")

peakLinks$sigDisease <- peakLinks$CorrelationDisease >= corCutOff & peakLinks$FDRDisease <= fdrCutOff 
peakLinks$sigHealthy <- peakLinks$CorrelationHealthy >= corCutOff & peakLinks$FDRHealthy <= fdrCutOff 


linksSig <- peakLinks[which(peakLinks$sigDisease | peakLinks$sigHealthy)]
linksDisease <- peakLinks[which(peakLinks$sigDisease & !peakLinks$sigHealthy)]
linksHealthy <- peakLinks[which(!peakLinks$sigDisease & peakLinks$sigHealthy)]
linksShared <- peakLinks[which(peakLinks$sigDisease & peakLinks$sigHealthy)]

outMatch <- list(
  seA = seA[unique(mcols(linksSig)$peakName),], 
  seB = seB[unique(mcols(linksSig)$gene_name),], 
  linksDisease = linksDisease,
  linksHealthy = linksHealthy,
  linksShared = linksShared,
  linksSig = linksSig,
  linksAll = peakLinks
  
)

saveRDS(outMatch,file.path(args$outdir,"Save_DiseaseHealthy_p2gLinks.rds"))

