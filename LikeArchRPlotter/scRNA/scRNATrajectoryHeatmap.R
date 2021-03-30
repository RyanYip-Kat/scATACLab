	
getGroupMatrix<-function (seurat=NULL, features = NULL, groupList = NULL,
    threads = 4,  verbose = TRUE,slot="data",splitBy=NULL,
    asSparse = FALSE,useFM=FALSE)
{
    require(ArchR)
    require(stringr)
    
    Donors <- unlist(lapply(Cells(seurat),function(x)return(str_split(x,"-")[[1]][2])))
    seurat$Donors <- Donors
    meta <- seurat@meta.data
    Donors <- unique(Donors)
    #seqnames <- unique(featureDF$seqnames)
    if(!useFM){
	    if(is.null(features) & !is.null(splitBy)){
		    features <- rownames(seurat)
		    Groups <- table(meta[[splitBy]])
		    features <- lapply(names(Groups),function(i)return(sample(features,size=Groups[[i]],replace=TRUE)))
	    }
    }else{
	    if(is.null(splitBy)){stop("groupBy must not be NULL!")}
	    require(future)
	    message("Find Markers ...")
	    plan("multiprocess", workers = threads)
	    Idents(seurat)=meta[[splitBy]]
	    markers <- FindAllMarkers(object=seurat,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
		       only.pos=FALSE,
                       slot=slot)
	    Groups <- unique(markers$cluster)
	    features <- lapply(Groups,function(i){
				       f <- markers$gene[markers$cluster%in%i]
				       return(f)})
    }
    #rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
    allMat <- GetAssayData(seurat,slot=slot)
    cellNames <- unlist(groupList, use.names = FALSE)
    allCellsList <- lapply(seq_along(Donors), function(x) {
        allCells <- rownames(meta)[meta$Donors%in%x]
        allCells <- allCells[allCells %in% cellNames]
        if (length(allCells) != 0) {
            allCells
        }
        else {
            NULL
        }
    })
    mat <- ArchR:::.safelapply(seq_along(features), function(x) {
        featureDFx <- features[[x]]
        matChr <- matrix(0, nrow =length(featureDFx), ncol = length(groupList))
        colnames(matChr) <- names(groupList)
        rownames(matChr) <- featureDFx
        for (y in seq_along(Donors)) {
            allCells <- allCellsList[[y]]
            if (!is.null(allCells)) {
	        maty <- allMat[featureDFx,allCells,drop=FALSE]
                for (z in seq_along(groupList)) {
                  cellsGroupz <- groupList[[z]]
                  idx <- BiocGenerics::which(colnames(maty) %in%
                    cellsGroupz)
                  if (length(idx) > 0) {
                    matChr[, z] <- matChr[, z] + Matrix::rowSums(maty[,
                      idx, drop = FALSE])
                  }
                }
                rm(maty)
            }
            if (y%%20 == 0 | y%%length(Donors) == 0) {
                gc()
            }
        }
        if (asSparse) {
            matChr <- as(matChr, "dgCMatrix")
        }
        cat(sprintf("Finished Group Matrix %s of %s\n",x, length(features)))
        matChr
    }, threads = threads) %>% Reduce("rbind", .)
    mat <- mat[unique(unlist(features)), , drop = FALSE]
    message("Successfully Created Group Matrix")
    gc()
    return(mat)
}

getGroupList <- function(
  seurat = NULL,
  name = "Trajectory",
  groupEvery = 1
  ){
  
  meta <-  seurat@meta.data
  trajectory <- meta[,name,drop=FALSE]
  trajectory <- trajectory[!is.na(trajectory[,1]),,drop=FALSE]
  breaks <- seq(0, 100, groupEvery)
  if(!all(is.numeric(trajectory[,1]))){
    stop("Trajectory must be a numeric. Did you add the trajectory with addTrajectory?")
  }
  if(!all(trajectory[,1] >= 0 & trajectory[,1] <= 100)){
    stop("Trajectory values must be between 0 and 100. Did you add the trajectory with addTrajectory?")
  }

  groupList <- lapply(seq_along(breaks), function(x){
      if(x == 1){
          NULL
      }else{
          rownames(trajectory)[which(trajectory[,1] > breaks[x - 1] & trajectory[,1] <= breaks[x])]
      }
  })[-1]
  names(groupList) <- paste0("T.", breaks[-length(breaks)], "_", breaks[-1])
  return(groupList)
}


CustomGetTrajectory <- function(
  seurat = NULL,
  splitBy=NULL,
  useFM=FALSE, # use FindAllMarkers
  name = "Trajectory",
  groupEvery = 1,
  log2Norm = TRUE,
  slot="data",
  features=NULL,
  scaleTo = 10000,
  smoothWindow = 11,
  threads = 16
  ){    
	message("Creating Trajectory Group List ..")
	groupList <- getGroupList(seurat,name=name,groupEvery=groupEvery)
	message("Creating Trajectory Group Matrix..")
	groupMat <- getGroupMatrix(seurat,
				   features=features,
				   useFM=useFM,
				   slot=slot,
				   groupList=groupList,
				   splitBy=splitBy,
				   threads=threads)
	if(!is.null(scaleTo)){
		if(any(groupMat < 0)){
		       	message("Some values are below 0, this could be a DeviationsMatrix in which scaleTo should be set = NULL.\nContinuing without depth normalization!")
		}else{
			groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
	       	}
	}
	
	if(log2Norm){
		if(any(groupMat < 0)){
			message("Some values are below 0, this could be a DeviationsMatrix in which log2Norm should be set = FALSE.\nContinuing without log2 normalization!")
		}else{
			groupMat <- log2(groupMat + 1)
		}
	}
	
	if(!is.null(smoothWindow)){
		message("Smoothing...")
                smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) .centerRollMean(x, k = smoothWindow))))
                colnames(smoothGroupMat) <- paste0(colnames(groupMat))
                colnames(groupMat) <- paste0(colnames(groupMat))

                #Create SE
		seTrajectory <- SummarizedExperiment(
						     assays = SimpleList(
									 smoothMat = as.matrix(smoothGroupMat),
                                                                         mat = as.matrix(groupMat)))
	}else{
		colnames(groupMat) <- paste0(colnames(groupMat))
		#Create SE
                seTrajectory <- SummarizedExperiment(assays = SimpleList(mat = as.matrix(groupMat)))
	}
	
	metadata(seTrajectory)$Params <- list(matrixClass="matrix",
					      scaleTo = scaleTo,
					      log2Norm = log2Norm,
					      smoothWindow = smoothWindow,
					      date = Sys.Date())
	seTrajectory

}


.centerRollMean <- function (v = NULL, k = NULL)
{
    o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
    if (k%%2 == 0) {
        o2 <- c(rep(o1[k], floor(k/2) - 1), o1[-seq_len(k - 1)],
            rep(o1[length(o1)], floor(k/2)))
    }
    else if (k%%2 == 1) {
        o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k - 1)],
            rep(o1[length(o1)], floor(k/2)))
    }
    else {
        stop("Error!")
    }
    o2
}


######################  main function
#splitBy="label_fine"
#se=CustomGetTrajectory(seurat,splitBy=splitBy)
#plotTrajectoryHeatmap(se)
