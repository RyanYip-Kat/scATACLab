LatentaddUMAP <- function(
  ArchRProj = NULL,
  latent=NULL, 
  name = "UMAP",
  nNeighbors = 40,
  minDist = 0.4,
  metric = "cosine",
  dimsToUse = NULL,
  sampleCells = NULL,
  outlierQuantile = 0.9,
  saveModel = TRUE,
  verbose = TRUE,
  seed = 1,
  force = FALSE,
  threads = 1,
  ...
  ){

  #############################################################################################
  # Default Parameters for Input Embeddings!
  #############################################################################################
  require(ArchR)
  #Merge Parameters
  embeddingParams <- list(...)
  nDims <- ncol(latent)
  stopifnot(nCells(ArchRProj)==nrow(latent))
  if(!is.null(dimsToUse)){
	  if(dimsToUse <= nDims){
		  X <- latent[,1:dimsToUse]
	  }else{
		  X <- latent
	  }
  }else{
	  X <- latent
  }
  embeddingParams$X <- X
  embeddingParams$n_neighbors <- nNeighbors
  embeddingParams$min_dist <- minDist
  embeddingParams$verbose <- verbose
  embeddingParams$metric <- metric

  estimateUMAP <- FALSE
  projectDF <- DataFrame(row.names = rownames(embeddingParams$X), projected = rep(0, nrow(embeddingParams$X))) #Projection ID
  if(!is.null(sampleCells)){
    if(sampleCells < nrow(embeddingParams$X)){
      message("Creating an Estimated UMAP by sub-sampling cells N = ", sampleCells, "!")
      saveModel <- TRUE
      idx <- sample(seq_len(nrow(embeddingParams$X)), sampleCells)
      cellNames <- rownames(embeddingParams$X)
      saveX <- embeddingParams$X[-idx, , drop = FALSE]
      embeddingParams$X <- embeddingParams$X[idx, , drop = FALSE]
      estimateUMAP <- TRUE
      projectDF[idx, 1] <- 1
    }
  }

  if(saveModel){
    embeddingParams$ret_nn <- TRUE
    embeddingParams$ret_model <- TRUE 
  }else{
    embeddingParams$ret_nn <- FALSE
    embeddingParams$ret_model <- FALSE      
  }

  #############################################################################################
  # Run Embedding
  #############################################################################################
  #Seed
  set.seed(seed)
  uwot_umap <- do.call(uwot::umap, embeddingParams)

  if(estimateUMAP){
    uwot_umap2 <- uwot::umap_transform(X = saveX, model = uwot_umap, n_threads = as.integer(threads), verbose = verbose)
    #We should check the distances
    knnRef <- as.vector(nabor::knn(data = uwot_umap[[1]], query = uwot_umap[[1]], k = 2)$nn.dists[,-1])
    knnProj <- as.vector(nabor::knn(data = uwot_umap[[1]], query = uwot_umap2, k = 1)$nn.dists)
    idxExclude <- which(knnProj >= quantile(knnRef, outlierQuantile))
    uwot_umap2[idxExclude, ] <- NA
  }

  #############################################################################################
  # Add Embedding to Project
  #############################################################################################
  nc <- ncol(embeddingParams$X)
  nr <- nrow(embeddingParams$X)

  if(saveModel){
    dir.create(file.path(getOutputDirectory(ArchRProj), "Embeddings"), showWarnings = FALSE)
    modelFile <- .tempfile(
      pattern = "Save-Uwot-UMAP-Params-Latent", 
      tmpdir = file.path(getOutputDirectory(ArchRProj), "Embeddings"),
      fileext = ".tar",
      addDOC = TRUE
    )
    #file.path(getOutputDirectory(ArchRProj), "Embeddings", paste0("Save-Uwot-UMAP-Params-",reducedDims,"-",.randomStr(),".tar"))
    saveModelTmp <- .saveUWOT(uwot_umap, modelFile)
    if(!file.exists(modelFile)){
      warning("Model was not saved properly, continuing without saving model!")
      modelFile <- NA
    }
    dfEmbedding <- data.frame(uwot_umap[[1]])
    colnames(dfEmbedding) <- paste0("Latent","#UMAP_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(embeddingParams$X)
    embeddingParams$X <- NULL

    if(estimateUMAP){
      dfEmbedding2 <- data.frame(uwot_umap2)
      colnames(dfEmbedding2) <- paste0("Latent","#UMAP_Dimension_",seq_len(ncol(dfEmbedding2)))
      rownames(dfEmbedding2) <- rownames(saveX)
      rm(uwot_umap2)
      dfEmbedding <- rbind(dfEmbedding, dfEmbedding2)
      dfEmbedding <- dfEmbedding[cellNames,,drop=FALSE]
    }

    ArchRProj@embeddings[[name]] <- SimpleList(
      df = dfEmbedding, 
      params = c(
        embeddingParams,
        dimsToUse = dimsToUse,
        #scaleDims = scaleDims,
        #corCutOff = corCutOff,
        nr=nr,
        nc=nc,
        uwotModel = modelFile,
        estimateUMAP = estimateUMAP,
        projectID = projectDF
      )
    )
  }else{
    dfEmbedding <- data.frame(uwot_umap)    
    colnames(dfEmbedding) <- paste0("Latent","#UMAP_Dimension_",seq_len(ncol(dfEmbedding)))
    rownames(dfEmbedding) <- rownames(embeddingParams$X)
    embeddingParams$X <- NULL
    ArchRProj@embeddings[[name]] <- SimpleList(
      df = dfEmbedding, 
      params = c(
        embeddingParams,
        dimsToUse = dimsToUse,
        nr=nr,
        nc=nc,
        uwotModel = NA,
        estimateUMAP = estimateUMAP,
        projectID = projectDF
      )
    )
  }   

  return(ArchRProj)

}

#New Save UWOT
.saveUWOT <- function(model, file){
  tryCatch({
    uwot::save_uwot(model = model, file = file, verbose = TRUE)
  }, error = function(e){
    .saveUWOT_Deprecated(model = model, file = file) #backwards to previous version
  })
}

#save_uwot does not work because tarring doesnt work for some reason on Stanford's compute server
#Adapted from save_uwot
.saveUWOT_Deprecated <- function(model, file){
  file <- file.path(normalizePath(dirname(file)), basename(file))
  wd <- getwd()
  mod_dir <- tempfile(pattern = "dir")
  dir.create(mod_dir)
  uwot_dir <- file.path(mod_dir, "uwot")
  dir.create(uwot_dir)
  model_tmpfname <- file.path(uwot_dir, "model")
  .safeSaveRDS(model, file = model_tmpfname)
  metrics <- names(model$metric)
  n_metrics <- length(metrics)
  for (i in seq_len(n_metrics)) {
      nn_tmpfname <- file.path(uwot_dir, paste0("nn", i))
      if (n_metrics == 1) {
          model$nn_index$save(nn_tmpfname)
          model$nn_index$unload()
          model$nn_index$load(nn_tmpfname)
      }
      else {
          model$nn_index[[i]]$save(nn_tmpfname)
          model$nn_index[[i]]$unload()
          model$nn_index[[i]]$load(nn_tmpfname)
      }
  }
  setwd(mod_dir)
  system2("tar", "-cvf uwot.tar uwot", stdout = NULL, stderr = NULL)
  o <- .fileRename("uwot.tar", file)
  setwd(wd)
  if (file.exists(mod_dir)) {
      unlink(mod_dir, recursive = TRUE)
  }
  return(o)
}

#New Save UWOT
.loadUWOT <- function(file){
  tryCatch({
    uwot::load_uwot(file = file, verbose = TRUE)
  }, error = function(e){
    .loadUWOT_Deprecated(file = file, nDim = nDim) #backwards to previous version
  })
}

#Adapted from load_uwot
.loadUWOT_Deprecated <- function(file, nDim = NULL){
    model <- NULL
    tryCatch({
        mod_dir <- tempfile(pattern = "dir")
        dir.create(mod_dir)
        utils::untar(file, exdir = mod_dir)
        model_fname <- file.path(mod_dir, "uwot/model")
        if (!file.exists(model_fname)) {
            stop("Can't find model in ", file)
        }
        model <- readRDS(file = model_fname)
        metrics <- names(model$metric)
        n_metrics <- length(metrics)
        for (i in seq_len(n_metrics)){
            nn_fname <- file.path(mod_dir, paste0("uwot/nn", i))
            if (!file.exists(nn_fname)) {
                stop("Can't find nearest neighbor index ", nn_fname, " in ", file)
            }
            metric <- metrics[[i]]
            if(length(model$metric[[i]]) == 0){
              if(!is.null(nDim)){
                nDim2 <- nDim
              }else{
                nDim2 <- length(model$metric[[i]])
              }
            }
            if(!is.null(nDim)){
              nDim2 <- nDim
            }
            ann <- uwot:::create_ann(metric, ndim = nDim2)
            ann$load(nn_fname)
            if (n_metrics == 1) {
                model$nn_index <- ann
            }else{
                model$nn_index[[i]] <- ann
            }
        }
    }, finally = {
        if (file.exists(mod_dir)) {
            unlink(mod_dir, recursive = TRUE)
        }
    })
    model 
}

.tempfile<-function (pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE)
{
    dir.create(tmpdir, showWarnings = FALSE)
    if (addDOC) {
        doc <- paste0("-Date-", Sys.Date(), "_Time-", gsub(":",
            "-", stringr::str_split(Sys.time(), pattern = " ",
                simplify = TRUE)[1, 2]))
    }
    else {
        doc <- ""
    }
    tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir,
        fileext = paste0(doc, fileext))
}

###############################
library(argparse)
library(ArchR)
parser <- ArgumentParser(description='ChromVAR handle motif from ArchR')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="ArchR Project Path")


parser$add_argument("--latent",
                    type="character",
                    default=NULL,
                    help="Embed DF csv from SCALE.py")

parser$add_argument("--name",
                    type="character",
                    default="SCALE_UMAP",
                    help="DF slot name")

#parser$add_argument("--outdir",
#                    type="character",
#                    default="output")
args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

#outDir=args$outdir
#makedir(outDir)

message("INFO : Loading dataset")
proj = loadArchRProject(args$project)
latent = read.csv(args$latent,row.names=1,stringsAsFactors=F,header=T)

message("INFO : Add Latent UMAP")
proj = LatentaddUMAP(proj,
		     latent=latent,
		     name=args$name,
		     dimsToUse=NULL)

message("INFO :  Save ")
saveRDS(proj,file.path(getOutputDirectory(proj),"Save-ArchR-Project.rds"))
message("INFO : Done")

