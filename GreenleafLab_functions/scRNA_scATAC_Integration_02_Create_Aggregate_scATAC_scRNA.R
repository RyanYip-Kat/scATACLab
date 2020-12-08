#Clustering and scATAC-seq UMAP for Hematopoiesis data
#06/02/19
#Cite Granja*, Klemm*, Mcginnis* et al. 
#A single cell framework for multi-omic analysis of disease identifies 
#malignant regulatory signatures in mixed phenotype acute leukemia (2019)
#Created by Jeffrey Granja
library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(argparse)
library(yaml)
library(Rcpp)
set.seed(1)

#functions="/home/ye/Work/R/scATAC/ArchR/GreenleafLab_functions/functions.R"
#source(functions)

####################################################
#Functions
####################################################

getGeneGTF <- function(file){
  #Import
  message("Reading in GTF...")
  importGTF <- rtracklayer::import(file)
  #Exon Info
  message("Computing Effective Exon Lengths...")
  exonGTF <- importGTF[importGTF$type=="exon",]
  exonList <- GenomicRanges::reduce(split(exonGTF, mcols(exonGTF)$gene_id))
  exonReduced <- unlist(exonList, use.names=TRUE)
  mcols(exonReduced)$gene_id <- names(exonReduced)
  mcols(exonReduced)$widths <- width(exonReduced)
  exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
  exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>% 
    unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
  #Gene Info
  message("Constructing gene GTF...")
  geneGTF1 <- importGTF[importGTF$type=="gene",]
  geneGTF2 <- GRanges(
      seqnames=paste0("chr",seqnames(geneGTF1)),
      ranges=ranges(geneGTF1),
      strand=strand(geneGTF1),
      gene_name=geneGTF1$gene_name,
      gene_id=geneGTF1$gene_id
    ) %>% keepFilteredChromosomes %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
  mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
  return(geneGTF2)
}

keepFilteredChromosomes <- function(x, remove = c("chrM"), underscore = TRUE, standard = TRUE, pruning.mode="coarse"){
	#first we remove all non standard chromosomes
	if(standard){
		x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = pruning.mode)
        }
        #Then check for underscores or specified remove
        seqNames <- seqlevels(x)
        chrRemove <- c()
        #first we remove all chr with an underscore
        if(underscore){
		chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
        }
        #next we remove all chr specified in remove
        chrRemove <- c(chrRemove, which(seqNames %in% remove))
        if(length(chrRemove) > 0){
		chrKeep <- seqNames[-chrRemove]
        }else{
		chrKeep <- seqNames
        }
        #this function restores seqlevels
        seqlevels(x, pruning.mode=pruning.mode) <- chrKeep
        return(x)
}
featureToGR <- function(feature){
  featureSplit <- stringr::str_split(paste0(feature), pattern = "_", n = 3, simplify = TRUE)
  GRanges(
    seqnames = featureSplit[,1],
    ranges = IRanges(
      start = as.integer(featureSplit[,2]),
      end = as.integer(featureSplit[,3]))
  )
}

grToFeature <- function(gr){
  paste(seqnames(gr),start(gr),end(gr),sep="_")
}

####################################################
#Input Data
####################################################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="./output")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

#Input Files( control and treatment combined object)
seATAC_file <- "./CEpi_scATAC/0h-3h/PeakMatrix_Summarized-Experiment.rds"
ciceroKNN_file <- "./CEpi_scATAC_Cicero/0h-3h/save-cicero-KNN-Groupings-cds.rds"
ciceroATAC_file <- "./CEpi_scATAC_Cicero/0h-3h/save-cicero-aggregated-accessibility-cds.rds"
seRNA_file <- "./CEpi_scRNA/0h-3h/scRNA_SummarizedExperiment_Expriment.rds"
CCA_file <- "./combined_scRNA_scATAC/0h-3h/Save-combined-KNN-UMAP.rds"
gtf_file <- "/home/ye/Data/10X/VDJ/ref/refdata-cellranger-mm10-3.0.0/genes/genes.gtf"
####################################################
#Get Clusters information for each KNN Group top group/exp wins!
scATAC <- readRDS(seATAC_file)
KNN <- data.frame(readRDS(ciceroKNN_file), stringsAsFactors=FALSE)
KNN <- apply(KNN,2,paste0)
KNNClusters <- apply(KNN, 2, function(x) colData(scATAC)[x,"Clusters"])
KNNGroups <- apply(KNN, 2, function(x) colData(scATAC)[x,"Group"])  # Groups' level ==2
KNN_Highest_Cluster <- lapply(seq_len(nrow(KNN)), function(x) names(sort(table(KNNClusters[x,]), decreasing=TRUE))[1]) %>% unlist
KNN_Highest_Experiment <- lapply(seq_len(nrow(KNN)), function(x) names(sort(table(KNNGroups[x,]), decreasing=TRUE))[1]) %>% unlist

####################################################
# scATAC-seq Merging from Cicero
####################################################
ciceroObj <- readRDS(ciceroATAC_file)
se <- SummarizedExperiment(
  assays = SimpleList(counts = assayData(ciceroObj)$exprs),
  rowRanges = featureToGR(featureData(ciceroObj)[[1]]),
  colData = DataFrame(
    row.names = colnames(assayData(ciceroObj)$exprs), 
    clustATAC = KNN_Highest_Cluster, 
    groupATAC = KNN_Highest_Experiment
  )
  )
metadata(se)$knn <- KNN
metadata(se)$knnClust <- KNNClusters
metadata(se)$knnGroup <- KNNGroups
rownames(se) <- grToFeature(rowRanges(se))
saveRDS(se, file.path(args$outdir,"Save-scATAC-Merged-KNN-SVD.rds"))  # control and treatment combined object

####################################################
# scRNA-seq Merging from Cicero
####################################################
scRNA <- readRDS(seRNA_file)
CCA <- readRDS(CCA_file)
CCA <- CCA$matchedCells
CCA <- CCA[CCA$corCCA >= 0.4,] #Filter low correlation to CCA R > 0.4

#Get scRNA Matrix
scMat <- assay(scRNA)

#Create aggregated Matched RNA Matrix exclude non-mappings
matRNA <- matrix(NA, ncol = nrow(KNN), nrow = nrow(scMat))
for(i in seq_len(nrow(KNN))){
  if(i%%100==0) print(i)
  knnIdx <- paste0(t(KNN[i,]))
  rnaIdx <- CCA$y[match(knnIdx, CCA$x)]
  rnaIdx <- na.omit(rnaIdx)
  matRNA[, i] <- Matrix::rowSums(scMat[,rnaIdx,drop=FALSE])
}
colnames(matRNA) <- colnames(se)

#Create aggregated summarized experiment
seRNA <- SummarizedExperiment(
    assays = SimpleList(counts = matRNA)
  )
rownames(seRNA) <- rownames(scMat)
gtf <- getGeneGTF(gtf_file)
gtfMatch <- gtf[!is.na(match(rownames(seRNA),gtf$gene_name))]
gtfMatch<-gtfMatch[gtfMatch$gene_name%in%rownames(seRNA)]

gtfMatch<-gtfMatch[!duplicated(gtfMatch$gene_name)]
seRNA<-seRNA[gtfMatch$gene_name,]
names(gtfMatch) <- rownames(seRNA)
rowRanges(seRNA) <- gtfMatch

#Use ATAC cluster info
colData(seRNA) <- DataFrame(row.names = colnames(seRNA), 
  clustATAC = KNN_Highest_Cluster, 
  groupATAC = KNN_Highest_Experiment
  )

saveRDS(seRNA, file.path(args$outdir,"Save-scRNA-Merged-KNN-SVD.rds"))

