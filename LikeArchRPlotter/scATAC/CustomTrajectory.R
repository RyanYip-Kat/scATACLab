.getQuantiles <- function (v = NULL, len = length(v))
{
    if (length(v) < len) {
        v2 <- rep(0, len)
        v2[seq_along(v)] <- v
    }
    else {
        v2 <- v
    }
    p <- trunc(rank(v2))/length(v2)
    if (length(v) < len) {
        p <- p[seq_along(v)]
    }
    return(p)
}


getSeuratTrajectory <- function(
  object = NULL,
  name = "Trajectory",
  trajectory = NULL, 
  groupBy = "Clusters",
  #reducedDims = "IterativeLSI",
  embedding = NULL,
  preFilterQuantile = 0.9, 
  postFilterQuantile = 0.9,
  useAll = FALSE, 
  dof = 250,
  spar = 1,
  force = FALSE,
  seed = 1
  ){
  
  stopifnot(class(object)=="Seurat")
  if(!is.null(seed)) set.seed(seed)
  meta.data <- object@meta.data
  stopifnot(groupBy%in%colnames(meta.data))
  groupDF <- meta.data[,groupBy,drop=FALSE] 
  groupDF <- groupDF[groupDF[,1] %in% trajectory,,drop=FALSE]

  if(sum(unique(groupDF[,1]) %in% trajectory)==0){
    stop("trajectory does not span any groups in groupBy! Are you sure your input is correct?")
  }

  if(sum(unique(groupDF[,1]) %in% trajectory) < 3){
    if(!force){
      stop("trajectory must span at least 3 groups in groupBy!")
    }
  }

  if(is.null(embedding)){
    mat <- Embeddings(object,reduction ="pca")
  }else{
    mat <- Embeddings(object, reduction = embedding)
  }
  mat <- mat[rownames(groupDF),,drop = FALSE]

  ######################################################
  #Filter Outliers
  ######################################################
  filterObj <- lapply(seq_along(trajectory), function(x){
      
      #Subset
      groupsx <- rownames(groupDF)[groupDF[,1]==trajectory[x]]
      matx <- mat[groupsx,,drop = FALSE]

      #Filter Distance
      matMeanx <- colMeans(matx)
      diffx <- sqrt(colSums((t(matx) - matMeanx)^2))
      idxKeep <- which(diffx <= quantile(diffx, preFilterQuantile))
      
      #Filter
      list(mat = matx[idxKeep,,drop=FALSE], groups = groupsx[idxKeep])

  })

  matFilter <- lapply(seq_along(filterObj), function(x) filterObj[[x]]$mat) %>% Reduce("rbind", .)
  groupsFilter <- groupDF[lapply(seq_along(filterObj), function(x) filterObj[[x]]$groups) %>% Reduce("c", .),,drop=FALSE]

  ######################################################
  #Now Initial Alignment
  ######################################################
  message("INFO : Initial Alignment Before Spline Fit")
  initialTime <- lapply(seq_along(trajectory), function(x){
      
      groupsx <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x]]
      matx <- matFilter[groupsx,,drop = FALSE]
      
      #Get Differences
      if(x != length(trajectory)){
          groupsxp1 <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x + 1]]
          meanx <- colMeans(matFilter[groupsxp1,,drop = FALSE])
          diffx <- sqrt(colSums((t(matx) - meanx)^2))
          timex <- (1 - .getQuantiles(diffx)) + x
      }else{
          groupsxm1 <- rownames(groupsFilter)[groupsFilter[,1] == trajectory[x - 1]]
          meanx <- colMeans(matFilter[groupsxm1,,drop = FALSE])
          diffx <- sqrt(colSums((t(matx) - meanx)^2))
          timex <- .getQuantiles(diffx) + x
      }
      
      timex

  }) %>% unlist

  ######################################################
  #Fit Cubic Splines
  ######################################################
  message("INFO :Spline Fit")
  matSpline <- lapply(seq_len(ncol(matFilter)), function(x){
    tryCatch({
      stats::smooth.spline(
          x = initialTime, 
          y = matFilter[names(initialTime), x], 
          df = dof, 
          spar = spar
      )[[2]]
    }, error = function(e){
      errorList <- list(
        it = x,
        x = initialTime, 
        y = matFilter[names(initialTime), x], 
        df = dof, 
        spar = spar
      )
      #.logError(e, fn = "smooth.spline", info = "", errorList = errorList, logFile = logFile)      
    })
  }) %>% Reduce("cbind",.) %>% data.frame()

  ######################################################
  # 1. KNN Fit vs Actual
  ######################################################
  message("INFO : KNN to Spline")
  knnObj <- nabor::knn(
      data =  matSpline,
      query = mat, 
      k = 3
  )

  #Estimate place along trajectory
  knnIdx <- knnObj[[1]]
  knnDist <- knnObj[[2]]
  knnDiff <- ifelse(knnIdx[,2] > knnIdx[,3], 1, -1)
  knnDistQ <- .getQuantiles(knnDist[,1])

  #Filter Outlier Cells to Trajectory for High Resolution
  idxKeep <- which(knnDist[,1] <= quantile(knnDist[,1], postFilterQuantile))
  dfTrajectory <- DataFrame(
      row.names = rownames(mat),
      Distance = knnDist[, 1],
      DistanceIdx = knnIdx[, 1] + knnDiff * knnDistQ
  )[idxKeep, , drop = FALSE]

  ######################################################
  # 2. Fit cells not in trajectory clusters
  ######################################################
  if(useAll){
    message("INFO : Aligning cells not in trajectory")
    
    if(is.null(embedding)){
      mat2 <- Embeddings(object, reduction = "pca")
    }else{
      mat2 <- Embeddings(object, reduction = embedding)
    }
    meta.data <- object@meta.data
    groupDF <- meta.data[,groupBy,] 
    groupDF <- groupDF[groupDF[,1] %ni% trajectory,,drop=FALSE]
    mat2 <- mat2[rownames(groupDF),,drop = FALSE]

    #Nearest Neighbors
    knnObj2 <- nabor::knn(
        data =  matSpline,
        query = mat2, 
        k = 3
    )

    #Estimate place along trajectory
    knnIdx2 <- knnObj2[[1]]
    knnDist2 <- knnObj2[[2]]
    knnDiff2 <- ifelse(knnIdx2[,2] > knnIdx2[,3], 1, -1)
    knnDistQ2 <- .getQuantiles(knnDist2[,1])

    #Keep Cells that are within the maximum distance of a cluster
    idxKeep <- which(knnDist2[,1] < max(dfTrajectory[,1]))
    dfTrajectory2 <- DataFrame(
        row.names = rownames(mat2),
        Distance = knnDist2[, 1],
        DistanceIdx = knnIdx2[, 1] + knnDiff2 * knnDistQ2
    )[idxKeep, , drop = FALSE]

    #Final Output
    dfTrajectory3 <- rbind(dfTrajectory, dfTrajectory2)
  }else{
    dfTrajectory3 <- dfTrajectory
  }
  
  dfTrajectory3$Trajectory <- 100 * .getQuantiles(dfTrajectory3[,2])
  tr <- as.data.frame(dfTrajectory3) 
  tr <- tr[,"Trajectory",drop=FALSE]
  object <- AddMetaData(object,tr,col.name=name)
  object
}



plotSeuratTrajectory <- function(
  object = NULL,
  embedding = "UMAP",
  trajectory = "Trajectory",
  name = "Trajectory",
  colorBy="cellcoldata",
  log2Norm = TRUE,
  pal = NULL,
  size = 0.2,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  quantHex = 0.5,
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 6,
  plotAs = NULL,
  smoothWindow = 5,
  addFit=NULL,
  addArrow=TRUE,
  ...
  ){

  require("ggplot2")
  require("ArchR")
  require("Seurat")
  if(is.null(quantCut)){
    quantCut <- c(0, 1)
  }


  ##############################
  # Plot Helpers
  ##############################
  .summarizeHex <- function(x = NULL){
    quantile(x, quantHex, na.rm = TRUE)
  }

  ##############################
  # Get Trajectory
  ##############################
  meta.data <- object@meta.data
  dfT <- meta.data[,trajectory,drop=FALSE]
  idxRemove <- which(is.na(dfT[,1]))

  ##############################
  # Get Embedding
  ##############################
  df <- as.data.frame(Embeddings(object, reduction = embedding))
  dfT <- cbind(df, dfT[rownames(df),])
  colnames(dfT) <- c("x", "y", "PseudoTime")

  #Parameters
  plotParams <- list(...)
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of  Trajectory")
  plotParams$baseSize <- baseSize
  if(is.null(addFit)){
	  plotParams$addFit="loess"
  }else{
	  plotParams$addFit=addFit
  }

  if(tolower(colorBy) == "coldata" | tolower(colorBy) == "cellcoldata"){
	  plotParams$color <- as.vector(meta.data[,name,drop=FALSE][rownames(df), 1])
          plotParams$discrete <- ArchR:::.isDiscrete(plotParams$color)
          plotParams$continuousSet <- "horizonExtra"
          plotParams$discreteSet <- "stallion"
          plotParams$title <- paste(plotParams$title, " colored by\ncolData : ", name)
          if(is.null(plotAs)){
		  plotAs <- "hexplot"
          }
  }else{
	  plotParams$continuousSet <- "solarExtra"
	  if(!is.null(log2Norm)){
		  #log2Norm <- TRUE
		  slot <- "data"
                  plotParams$continuousSet <- "horizonExtra"
	  }
	  if(is.null(log2Norm)){
		  #log2Norm <- FALSE
		  slot <- "counts"
          }
	  MatSc <- GetAssayData(object,slot)
	  color <- as.matrix(MatSc[name,,drop=FALSE])
	  plotParams$color <- color[, rownames(df), drop = FALSE]
	  plotParams$discrete <- FALSE
          plotParams$title <- sprintf("%s colored by\n%s : %s", plotParams$title, colorBy, name)
          if(is.null(plotAs)){
		  plotAs <- "hexplot"
	  }
          if(plotAs=="hexplot"){
		  plotParams$fun <- .summarizeHex
	  }

  }


  #Additional Params!
  plotParams$xlabel <- colnames(df)[1]
  plotParams$ylabel <- colnames(df)[2]

  if(!is.null(continuousSet)){
    plotParams$continuousSet <- continuousSet
  }
  if(!is.null(continuousSet)){
    plotParams$discreteSet <- discreteSet
  }
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize

  if(plotParams$discrete){
    plotParams$color <- paste0(plotParams$color)
  }

  if(!plotParams$discrete){

    plotParams$color <- as.vector(plotParams$color)

    if(name != trajectory){
      plotParams$color <- ArchR:::.quantileCut(plotParams$color, min(quantCut), max(quantCut))
    }

    if(!is.null(log2Norm)){
      if(log2Norm){
        plotParams$color <- log2(plotParams$color + 1)
        plotParams$colorTitle <- paste0("Log2Norm")
      }else{
        plotParams$colorTitle <- "Raw"
      }
    }

    plotParams$color[idxRemove] <- NA
    plotParams$pal <- paletteContinuous(set = plotParams$continuousSet)

    if(tolower(plotAs) == "hex" | tolower(plotAs) == "hexplot"){
      plotParams$addPoints <- TRUE
      if(is.null(plotParams$bins)){
        plotParams$bins <- 150
      }
      message("Plotting")
      out <- do.call(ggHex, plotParams)
    }else{
      message("Plotting")
      out <- do.call(ggPoint, plotParams)
    }

  }else{
    message("Plotting")
    out <- do.call(ggPoint, plotParams)
  }

  if(!keepAxis){
    out <- out + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  message("Plotting Trajectory")

  #Prep Trajectory Vector
  dfT$value <- plotParams$color
  dfT <- dfT[order(dfT$PseudoTime), ]
  dfT <- dfT[!is.na(dfT$PseudoTime), ]

  #Plot Pseudo-Time
  out2 <- ggPoint(
    x = dfT$PseudoTime,
    y = dfT$value,
    color = dfT$PseudoTime,
    discrete = FALSE,
    xlabel = "PseudoTime",
    ylabel = name,
    pal = plotParams$pal,
    ratioYX = 0.5,
    rastr = TRUE
  ) + geom_smooth(color = "black")

  attr(out2, "ratioYX") <- 0.5

  if(addArrow){

    dfArrow <-  split(dfT, floor(dfT$PseudoTime / 1.01)) %>%
      lapply(colMeans) %>% Reduce("rbind",.) %>% data.frame
    dfArrow$x <- ArchR:::.centerRollMean(dfArrow$x, smoothWindow)
    dfArrow$y <- ArchR:::.centerRollMean(dfArrow$y, smoothWindow)
    dfArrow <- rbind(dfArrow, dfT[nrow(dfT), ,drop = FALSE])

    out <- out + geom_path(
            data = dfArrow, aes(x, y, color=NULL), size= 1,
            arrow = arrow(type = "open", length = unit(0.1, "inches"))
          )
  }
  list(out, out2)

}


##################

#seurat=readRDS("/Path/to/seurat")
#trajectory=c("Memory BC","Naive BC","Plasma BC")
#groupBy="label_fine"
#seurat=getSeuratTrajectory(seurat,groupBy=groupBy,trajectory=trajectory,embedding="umap")
#p1=plotSeuratTrajectory(seurat,embedding="umap",name="PAX6",colorBy="matrix")
#ggsave("plot2.pdf",plot=p1[[2]])
#ggsave("plot1.pdf",plot=p1[[1]])

