library(argparse)
library(stringr)
library(Seurat)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(Rcpp)
set.seed(1)

functions="/home/ye/Work/R/scATAC/ArchR/GreenleafLab_functions/functions.R"
source(functions)

####################################################
#Functions
####################################################

#Nearest Neighbor differential
findNN <- function(query, reference, method = "euclidean"){
    findClosest <- function(x, m, method = "euclidean"){
        if(method=="euclidean"){
            which.min(sqrt(colSums((t(m) - x) * (t(m) - x))))
        }else if(method=="pearson"){
            which.max(cor(t(m),x,method = method)[,1])
        }else if(method=="spearman"){
            which.max(cor(t(m),x,method = method)[,1])
        }
    }
    pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
    mat <- data.frame(matrix(ncol = 4, nrow = nrow(query)))
    colnames(mat) <- c("x", "i", "y", "j")
    for(i in seq_len(nrow(query))){
        setTxtProgressBar(pb,round(i*100/nrow(query),0))
      j <- findClosest(query[i,], reference, method)
      mat[i,] <- c(x = rownames(query)[i], i = i, y = rownames(reference)[j], j = j)
    }
    return(mat)
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

####################################################
#Input Data
####################################################


parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="./output")


parser$add_argument("--input_RNA",
                    type="character",
                    default=NULL,
		    help="RNA se object")

parser$add_argument("--input_GS",
                    type="character",
                    default=NULL,
		    help="the gene score se of scATAC from ArchR")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("### Loading dataset")
se <- readRDS(args$input_RNA)
matRNA <- assay(se)

#Prep Gene Score Matrix from Summarized Experiment
seGS <- readRDS(args$input_GS)
matGS <- assay(seGS)


#Parameters
nCCA <- 20
nVarGenes <- 2500
selectMethod <- "all"

#Gene Universe
geneUniverse <- intersect(rownames(matGS),rownames(matRNA))

#Remove Mito RNA
#geneUniverse <- geneUniverse[geneUniverse %ni% grep("^MT", c(rownames(seGS),rownames(se)), value = TRUE)]

#Subset By Gene Universe
matRNA <- matRNA[geneUniverse, ,drop = FALSE]
matGS <- matGS[geneUniverse, ,drop = FALSE]

#Create RNA Seurat
objRNA <- CreateSeuratObject(matRNA, project = "RNA")
objRNA <- NormalizeData(object = objRNA)
objRNA <- ScaleData(object = objRNA)
objRNA <- FindVariableFeatures(object = objRNA,selection.method= "dispersion", nfeatures= as.integer(nVarGenes))
objRNA$protocol <- "RNA"

#Create GS Seurat
objGS <- CreateSeuratObject(matGS, project = "ATAC")
objGS <- NormalizeData(object = objGS)
objGS <- ScaleData(object = objGS)
objGS <- FindVariableFeatures(object = objGS,selection.method= "dispersion", nfeatures= as.integer(nVarGenes))
objGS$protocol <- "ATAC"

#Intersect Variable Genes
if(tolower(selectMethod) == "genescores"){
  varGenes <- VariableFeatures(objGS)
}else if(tolower(selectMethod) == "rna"){
  varGenes <- VariableFeatures(objRNA)
}else if(tolower(selectMethod) == "intersect"){
  varGenes <- intersect(VariableFeatures(objGS),VariableFeatures(objRNA))
}else if(tolower(selectMethod) == "all"){
  varGenes <- unique(c(VariableFeatures(objGS),VariableFeatures(objRNA)))
}

#Run combined Seurat v2.3.4
#print("### Run combined")
#combined <- RunCCA(object = objRNA, object2 = objGS, features = varGenes, num.cc = as.integer(nCCA))

#Variance Expectation Ration Seurat v2.3.4
#combined <- CalcVarExpRatio(object = combined, reduction.type = "pca", grouping.var = "protocol", dims.use = seq_len(as.integer(nCCA)))

#Filter Seurat v2.3.4
#combined <- SubsetData(object = combined, subset.name = "var.ratio.pca", accept.low = 0.5)

#Align Subspace Seurat v2.3.4
#combined <- AlignSubspace(object = combined, reduction.type = "cca", grouping.var = "protocol", dims.align = seq_len(as.integer(nCCA)))

#saveRDS(combined, file.path(args$outdir,"Save-CCA-Alignment-scATAC-scRNA.rds"))
anchors <- FindIntegrationAnchors(object.list = list(objRNA,objGS),anchor.features=varGenes,
				  dims = 1:30,reduction="cca")
combined <- IntegrateData(anchorset = anchors,features.to.integrate=varGenes)

DefaultAssay(combined) <- "integrated"
saveRDS(combined,file.path(args$outdir,"Alignment-scATAC-scRNA.rds"))
# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

#Get combined Matrix
#alignedcombined<-Embeddings(combined, reduction = "cca")
#alignedcombined <- GetCellEmbeddings(combined, reduction.type = "cca.aligned")
pca.combined<-Embeddings(combined, reduction = "pca")
#KNN Search
#Alternatively for speed FNN::getknnx(query, reference, k = 1)
#We just used a simple function
norm.data=GetAssayData(combined,"data")
hvg=VariableFeatures(combined)
matchedCells <- findNN(
  query = pca.combined[combined@meta.data$protocol=="ATAC",],
  reference = pca.combined[combined@meta.data$protocol=="RNA",],
  method = "euclidean")

matchedCells$corCCA <- rowCorCpp(
  match(matchedCells$x, colnames(combined)),
  match(matchedCells$y, colnames(combined)),
  pca.combined, pca.combined)

matchedCells$corVarRNA <- rowCorCpp(
  match(matchedCells$x, colnames(combined)),
  match(matchedCells$y, colnames(combined)),
  t(as.matrix(norm.data[hvg,])),
  t(as.matrix(norm.data[hvg,])))

matchx <- match(matchedCells$x, colnames(combined))
matchy <- match(matchedCells$y, colnames(combined))
mat <- as.matrix(norm.data[hvg,])

#-------------------------------------------------------
#UMAP
#-------------------------------------------------------
set.seed(1)
umap <- uwot::umap(
    pca.combined,
    n_neighbors = 50,
    min_dist = 0.5,
    metric = "euclidean",
    n_threads = 5,
    verbose = TRUE,
    ret_model = FALSE)

#Plot DF
plotDF <- data.frame(umap)
rownames(plotDF) <- rownames(pca.combined)
plotDF[rownames(combined@meta.data[rownames(plotDF),]),"protocol"] <- combined@meta.data[rownames(plotDF),]$protocol
plotDF <- plotDF[sample(seq_len(nrow(plotDF)), nrow(plotDF)),, drop = FALSE]

saveRDS(list(plotDF = plotDF, matchedCells = matchedCells),file.path(args$outdir, "Save-combined-KNN-UMAP.rds"))
