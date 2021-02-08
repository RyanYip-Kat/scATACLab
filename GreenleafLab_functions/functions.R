library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(yaml)
library(Rcpp)
library(tidyverse)
library(uwot)
library(edgeR)
library(Seurat)
library(matrixStats)
library(cicero)
library(m3addon)
library(monocle)
#Helper function for summing sparse matrix groups
binarizeMat <- function(mat){
    mat@x[mat@x > 0] <- 1
    mat
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


sparseMatTTest <- function(mat1, mat2, m0 = 0){
	#Get Population Values
	n1 <- ncol(mat1)
	n2 <- ncol(mat2)
	n <- n1 + n2
	#Sparse Row Means
	m1 <- Matrix::rowMeans(mat1, na.rm=TRUE)
	m2 <- Matrix::rowMeans(mat2, na.rm=TRUE)
	#Sparse Row Variances
	#v1 <- ArchRx:::computeSparseRowVariances(mat1@i + 1, mat1@x, m1, n1)
	#v2 <- ArchRx:::computeSparseRowVariances(mat2@i + 1, mat2@x, m2, n2)
	v1=sparseRowVariances(mat1)
	v2=sparseRowVariances(mat2)
	#Calculate T Statistic
	se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*v1 + (n2-1)*v2)/(n1+n2-2) )
    tstat <- (m1-m2-m0)/se
	#tstat <- sqrt((n1 * n2) / n) / sqrt((n1-1)/(n-2)*v1 + (n2-1)/(n-2)*v2)
	pvalue <- 2*pt(-abs(tstat), n - 2)
	fdr <- p.adjust(pvalue, method = "fdr")
	out <- data.frame(fdr = fdr, pval = pvalue, tstat = tstat, mean1 = m1, mean2 = m2, var1 = v1, var2 = v2, n1 = n1, n2 = n2)
	return(out)
}

groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
        else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

#LSI Adapted from fly-atac with information for re-projection analyses
calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){

    set.seed(1)

    #TF IDF LSI adapted from flyATAC
    if(binarize){
        message(paste0("Binarizing matrix..."))
        mat@x[mat@x > 0] <- 1
    }

    #Filter 0 Sum Peaks
    rowSm <- Matrix::rowSums(mat)
    if(!is.null(nFeatures)){
        message(paste0("Getting top ", nFeatures, " features..."))
        idx1 <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
        idx2 <- which(rowSm>0)
        idx <- intersect(idx1,idx2)
        mat <- mat[idx,,drop=FALSE]
    }else{
        idx <- which(rowSm>0)
        mat <- mat[idx,,drop=FALSE]
    }

    #Filter 0 Sum Cells
    colSm <- Matrix::colSums(mat)
    if(length(which(colSm==0))>0){
        message("Filtering Cells with 0 ColSums...")
        mat <- mat[,which(colSm>0),drop=FALSE]
    }

    #Calc RowSums and ColSums
    colSm <- Matrix::colSums(mat)
    rowSm <- Matrix::rowSums(mat)

    #Calc TF IDF
    message("Computing Term Frequency IDF...")
    freqs <- t(t(mat)/colSm)
    idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs

    #Calc SVD then LSI
    message("Computing SVD using irlba...")
    svd <- irlba::irlba(tfidf, nComponents, nComponents)
    svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))

    #Return Object
    out <- list(
        matSVD = matSVD,
        rowSm = rowSm,
        colSm = colSm,
        idx = idx,
        svd = svd,
        binarize = binarize,
        nComponents = nComponents,
        date = Sys.Date(),
        seed = 1)

    out

}

projectLSI <- function(mat, lsi){

    #Get Same Features
    mat <- mat[lsi$idx,]
    if(lsi$binarize){
        message(paste0("Binarizing matrix..."))
        mat@x[mat@x > 0] <- 1
    }

    #Calc TF IDF
    rowsToZero <- which(lsi$rowSm == 0)
    setToZero <- which((mat@i + 1) %in% rowsToZero)
    if(length(setToZero) > 0){
        mat@x[setToZero] <- 0
    }

    message("Computing Term Frequency IDF...")
    freqs <- t(t(mat)/Matrix::colSums(mat))
    idf   <- as(log(1 + length(lsi$colSm) / lsi$rowSm), "sparseVector")
    tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
    if(length(Matrix::which(is.na(tfidf),arr.ind=TRUE)) > 0){
        tfidf[Matrix::which(is.na(tfidf),arr.ind=TRUE)] <- 0 #weird Inf * 0
    }

    #Calc V
    V <- t(tfidf) %*% lsi$svd$u %*% diag(1/lsi$svd$d)

    #Calc SVD then LSI
    message("Computing SVD using irlba...")
    svdDiag <- matrix(0, nrow=lsi$nComponents, ncol=lsi$nComponents)
    diag(svdDiag) <- lsi$svd$d
    matSVD <- t(svdDiag %*% t(V))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))

    return(matSVD)

}

#Sparse Variances Rcpp
sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  using namespace std;
  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)

#Seurat SNN
seuratSNN <- function(matSVD, dims.use = 1:50, print.output = TRUE,resolution=0.8){
  set.seed(1)
  message("Making Seurat Object...")
  mat <- matrix(rnorm(nrow(matSVD) * 3, 1000), ncol = nrow(matSVD), nrow = 3)
  colnames(mat) <- rownames(matSVD)
  rownames(mat)=paste("M",1:nrow(mat))
  obj <- CreateSeuratObject(mat, project='scATAC', min.cells=0, min.features=0)
  obj[["pca"]] <- CreateDimReducObject(embeddings=matSVD,key="PCA_",assay = DefaultAssay(obj))
  obj <- FindNeighbors(obj,reduction="pca",dims=dims.use)
  obj <- FindClusters(object = obj, resolution=resolution,verbose=print.output)
  #clust <- obj@meta.data[,ncol(obj@meta.data)]
  #paste0("Cluster",match(clust, unique(clust)))
  clust=obj$seurat_clusters
  return(paste("Cluster",clust))
}

##########################
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

#############################
featureToGR <- function(feature){
    featureSplit <- stringr::str_split(paste0(feature), pattern = "_", n = 3, simplify = TRUE)
    gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
    return(gr)
}

###############################
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

#############################
getGeneGTF <- function(file){
    #Import
    message("Reading in GTF...")
    importGTF <- rtracklayer::import(file)
    #Exon Info
    message("Computing Effective Exon Lengths...")
    exonGTF <- importGTF[importGTF$type=="exon",]
    exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
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
# gtfFile <- "data/genes.gtf"
#genes <- getGeneGTF(gtfFile) %>% resize(1,"start") %>% resize(tssWindow * 2 + 1, "center")
###################################
#cmdPeaks <- sprintf(
#	    "macs2 callpeak -g %s --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all", 
#	    genome_size, 
#	    names(clusterResults)[j], 
#	    clusterBedj, 
#	    dirPeaks
#	  )


#########################
#Optimized LSI for scRNA-seq analysis
scRNA_optimizeLSI <- function(mat, scaleTo = 10000, priorCount = 3, pcsUse = 1:25,
    resolution = c(0.2, 0.4, 0.8), varFeatures = c(2500, 2500, 2500), seed = 1){

    set.seed(seed)
    stopifnot(length(resolution) > 1)

    #Initialize List
    lsiOut <- list()

    #Initial LSI uses variances that are across all single cells and will have larger batch relationships
    i <- 1
    message("Initial LSI...")
    matNorm <- t(t(mat)/Matrix::colSums(mat)) * scaleTo
    matNorm@x <- log2(matNorm@x + 1)
    idVarFeatures <- head(order(sparseRowVariances(matNorm),decreasing=TRUE), varFeatures[i])
    lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
    clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], print.output = FALSE)

    #Store
    lsiOut[[paste0("iter", i)]] <- list(
        lsiMat = lsiObj$matSVD,
        varFeatures = idVarFeatures,
        clusters = clusters
        )

    for(i in seq(2, length(varFeatures))){

       message(sprintf("Additional LSI %s...", i))

        #Run LSI
        clusterMat <- edgeR::cpm(groupSums(mat, clusters, sparse = TRUE), log=TRUE, prior.count = priorCount)
        idVarFeatures <- head(order(rowVars(clusterMat), decreasing=TRUE), varFeatures[i])
        lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
        clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], print.output = FALSE)

        if(i == length(varFeatures)){
            #Save All Information from LSI Attempt
            lsiOut[[paste0("iter", i)]] <- list(
                lsiObj = lsiObj,
                varFeatures = idVarFeatures,
                clusters = clusters,
                matNorm = matNorm
                )
        }else{
            lsiOut[[paste0("iter", i)]] <- list(
                lsiMat = lsiObj$matSVD,
                varFeatures = idVarFeatures,
                clusters = clusters
                )
        }

    }

    return(lsiOut)

}
