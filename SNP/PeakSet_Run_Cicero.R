#Running gwas chromVAR on single-cell summarized experiment using co-accessibility
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(ArchR)
library(argparse)
library(stringr)
library(cicero)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(yaml)
library(Rcpp)
set.seed(1)


parser <- ArgumentParser(description='Computing Gene Activity scores using Cicero and Co-Accessibility')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="the project path  of ArchR")


parser$add_argument("--outdir",
                    type="character",
                    default="./results")


parser$add_argument("--binarize",
                    action="store_true")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


grToFeature <- function(gr){
	peakinfo <- data.frame(
		row.names = paste(seqnames(gr),start(gr),end(gr),sep="_"),
		site_name = paste(seqnames(gr),start(gr),end(gr),sep="_"),
		chr = gsub("chr","",as.character(seqnames(gr))),
		bp1 = start(gr),
		bp2 = end(gr)
	)
	return(peakinfo)
}

featureToGR <- function(feature){
	featureSplit <- stringr::str_split(paste0(feature), pattern = "_", n = 3, simplify = TRUE)
	gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
	return(gr)
}

makeCDS <- function(se, binarize = TRUE){
	peakinfo <- grToFeature(se)
	mat <- assay(se)
	if(binarize){
		mat@x[which(mat@x > 0)] <- 1
	}
	cellinfo <- data.frame(colData(se))
	cellinfo$cells <- rownames(cellinfo)
	cds <-  suppressWarnings(newCellDataSet(mat,
	                          phenoData = methods::new("AnnotatedDataFrame", data = cellinfo),
	                          featureData = methods::new("AnnotatedDataFrame", data = peakinfo),
	                          expressionFamily=negbinomial.size(),
	                          lowerDetectionLimit=0))
	fData(cds)$chr <- as.character(fData(cds)$chr)
	fData(cds)$bp1 <- as.numeric(as.character(fData(cds)$bp1))
	fData(cds)$bp2 <- as.numeric(as.character(fData(cds)$bp2))
	cds <- cds[order(fData(cds)$chr, fData(cds)$bp1),]
	return(cds)
}


sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  using namespace std;
  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {

    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }
    if(max(idxX) > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }
    if(max(idxY) > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }
    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);

    const int nx = X.ncol();
    const int ny = Y.ncol();
    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }
    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }
    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }
    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }
    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);
    }
    return(cor);
  }'
)

getTxDbGenes <- function(txdb = NULL, orgdb = NULL, gr = NULL, ignore.strand = TRUE){

    if (is.null(genome)) {
        if (is.null(txdb) | is.null(orgdb)) {
            stop("If no provided genome then you need txdb and orgdb!")
        }
    }

    if (is.null(gr)) {
        genes <- GenomicFeatures::genes(txdb)
    }else {
        genes <- suppressWarnings(subsetByOverlaps(GenomicFeatures::genes(txdb), gr, ignore.strand = ignore.strand))
    }

    if (length(genes) > 1) {
        mcols(genes)$symbol <- suppressMessages(mapIds(orgdb,
            keys = mcols(genes)$gene_id, column = "SYMBOL", keytype = "ENTREZID",
            multiVals = "first"))
        genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
        names(genes) <- NULL
        out <- genes
    }else {
        out <- GRanges(seqnames(gr), ranges = IRanges(0, 0), gene_id = 0, symbol = "none")[-1]
    }

    return(out)

}

################################
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

message("INFO : Get Peak Set ")
projHeme=loadArchRProject(args$project)
obj=getMatrixFromProject(projHeme,useMatrix="PeakMatrix")
peaks=getPeakSet(projHeme)
rowRanges(obj)=peaks
# mcols(peaks)
rowData(obj)=cbind(rowData(obj),mcols(obj))
rgs=as.data.frame(ranges(peaks))
sqs=as.data.frame(seqnames(peaks))
peakset=cbind(sqs,rgs)
peak_set=paste(peakset$value,peakset$start,peakset$end,sep="_")
rownames(obj)=peak_set


#Read input
mdata <- colData(obj)
tssWindow <- 2500
flank <- 250*10^3
corCutOff <- 0.35
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#orgdb <- org.Hs.eg.db

#Reduced Dimensions
message("INFO : Reduced Dimensions ")
mat=as.data.frame(getEmbedding(projHeme,"UMAP"))
colnames(mat)=c("UMAP1","UMAP2")
colData(obj)$UMAP1=mat$UMAP1
colData(obj)$UMAP2=mat$UMAP2
dimred <- data.frame(
		row.names = colnames(obj),
		colData(obj)$UMAP1,
		colData(obj)$UMAP2
	)

#Get ChromSizes
message("INFO : Get ChromSizes")
chromSizes <- seqlengths(bsgenome)[paste0("chr",c(1:22,"X"))]
genome <- data.frame(names(chromSizes),chromSizes)
rownames(genome) <- NULL

message("INFO : Save scATAC-Summarized-Experiment")
saveRDS(obj,file.path(args$outdir,"scATAC-Summarized-Experiment.rds"))
#Get CDS
message("INFO : Get CDS")
obj <- makeCDS(obj, binarize = args$binarize)
obj <- detectGenes(obj)
obj <- estimateSizeFactors(obj)
ciceroObj <- make_cicero_cds(obj, k = 50, reduced_coordinates = dimred[colnames(obj),])

#Compute Correlations
message("INFO : Computing grouped correlations...")
gr <- featureToGR(featureData(ciceroObj)[[1]])
o <- suppressWarnings(as.matrix( findOverlaps(resize( resize(gr,1,"center"), 2*flank + 1, "center"), resize(gr,1,"center"), ignore.strand=TRUE) ))
o <- data.table::as.data.table(data.frame(i = matrixStats::rowMins(o), j = matrixStats::rowMaxs(o)))
o <- data.frame(o[!duplicated(o),])
o <- o[o[,1]!=o[,2],]
o$cor <- rowCorCpp(o[,1], o[,2], assayData(ciceroObj)$exprs, assayData(ciceroObj)$exprs)
connections <- data.frame(
	Peak1 = featureData(ciceroObj)[[1]][o[,1]],
	Peak2 = featureData(ciceroObj)[[1]][o[,2]],
	coaccess = o[,3]
	)

#Save Output
message("INFO : Save Output")
saveRDS(ciceroObj,file.path(args$outdir,"ciceroObj.rds"))
saveRDS(connections, file.path(args$outdir,"Peaks-Co-Accessibility.rds"))
