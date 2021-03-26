library(ArchR)
source("/home/ye/Work/R/scATAC/ArchR/plotter/ArchRHeatmap.R")
MyplotMarkers <- function(
  seMarker = NULL,
  name = NULL,
  useName="Pval",
  plotAs = "Volcano",
  scaleTo = 10^4
  ){

  #Evaluate AssayNames
  if(is.null(useName)){
	  stop("useName must be not NULL!")
  }

  cutOff=paste0(useName," <= 0.01 & abs(Log2FC) >= 0.5")
  assayList=assays(seMarker)
  assayNames <- names(SummarizedExperiment::assays(seMarker))
  for(an in assayNames){
    eval(parse(text=paste0(an, " <- ", "SummarizedExperiment::assays(seMarker)[['", an, "']]")))
  }
  passMat <- eval(parse(text=cutOff))
  for(an in assayNames){
    eval(parse(text=paste0("rm(",an,")")))
  }
  passMat[is.na(passMat)] <- FALSE

  if(is.null(name)){
    name <- colnames(seMarker)[1]
  }

  DF <- assays(seMarker[,name])[[useName]]
  DF <- as.vector(as.matrix(DF))
  DF[is.na(DF)] <- 1

  if(tolower(plotAs) == "volcanodiff"){
    Diff <- assays(seMarker[,name])$MeanDiff
    Diff <- as.vector(as.matrix(Diff))
    qDiff <- max(quantile(abs(Diff), probs = 0.999, na.rm=TRUE), 4) * 1.05
    color <- ifelse(passMat[, name], "Differential", "Not-Differential")
    color[color == "Differential"] <- ifelse(Diff[color == "Differential"] > 0, "Up-Regulated", "Down-Regulated")
  }else{
    LFC <- assays(seMarker[,name])$Log2FC
    LFC <- as.vector(as.matrix(LFC))
    qLFC <- max(quantile(abs(LFC), probs = 0.999, na.rm=TRUE), 4) * 1.05
    LM <- log2((assays(seMarker[,name])$Mean + assays(seMarker[,name])$MeanBGD)/2 + 1)
    LM <- as.vector(as.matrix(LM))
    color <- ifelse(passMat[, name], "Differential", "Not-Differential")
    color[color == "Differential"] <- ifelse(LFC[color == "Differential"] > 0, "Up-Regulated", "Down-Regulated")
  }

  pal <- c("Up-Regulated" = "firebrick3", "Not-Differential" = "lightgrey", "Down-Regulated" = "dodgerblue3")

  idx <- c(which(!passMat[, name]), which(passMat[, name]))
  title <- sprintf("Number of features = %s\nNumber Up-Regulated = %s (%s Percent)\nNumber Down-Regulated = %s (%s Percent)",
    nrow(seMarker), sum(color=="Up-Regulated"),
    round(100 * sum(color=="Up-Regulated") / nrow(seMarker), 2),
    sum(color=="Down-Regulated"),
    round(100 * sum(color=="Down-Regulated") / nrow(seMarker), 2)
    )

  if(tolower(plotAs) == "ma"){
    ggPoint(
      x = LM[idx],
      y = LFC[idx],
      color = color[idx],
      ylim = c(-qLFC, qLFC),
      size = 1,
      extend = 0,
      rastr = TRUE,
      labelMeans = FALSE,
      labelAsFactors = FALSE,
      pal = pal,
      xlabel = "Log2 Mean",
      ylabel = "Log2 Fold Change",
      title = title
    ) + geom_hline(yintercept = 0, lty = "dashed") +
    scale_y_continuous(breaks = seq(-100, 100, 2), limits = c(-qLFC, qLFC), expand = c(0,0))
  }else if(tolower(plotAs) == "volcano"){
    ggPoint(
      x = LFC[idx],
      y = -log10(DF[idx]),
      color = color[idx],
      xlim = c(-qLFC, qLFC),
      extend = 0,
      size = 1,
      rastr = TRUE,
      labelMeans = FALSE,
      labelAsFactors = FALSE,
      pal = pal,
      xlabel = "Log2 Fold Change",
      ylabel = paste0("-Log10 ",useName),
      title = title
    ) + geom_vline(xintercept = 0, lty = "dashed") +
    scale_x_continuous(breaks = seq(-100, 100, 2), limits = c(-qLFC, qLFC), expand = c(0,0))
  }else if(tolower(plotAs) == "volcanodiff"){
    ggPoint(
      x = Diff[idx],
      y = -log10(DF[idx]),
      color = color[idx],
      xlim = c(-qDiff, qDiff),
      extend = 0,
      size = 1,
      rastr = TRUE,
      labelMeans = FALSE,
      labelAsFactors = FALSE,
      pal = pal,
      xlabel = "Mean Difference",
      ylabel = "-Log10 useName",
      title = title
    ) + geom_vline(xintercept = 0, lty = "dashed")
  }else{
    stop("plotAs not recognized")
  }

}

computeKNN<-function (data = NULL, query = NULL, k = 50, includeSelf = FALSE,
    ...)
{
    if (is.null(query)) {
        query <- data
        searchSelf <- TRUE
    }
    else {
        searchSelf <- FALSE
    }
    require("nabor")
    if (searchSelf & !includeSelf) {
        knnIdx <- nabor::knn(data = data, query = query, k = k +
            1, ...)$nn.idx
        knnIdx <- knnIdx[, -1, drop = FALSE]
    }
    else {
        knnIdx <- nabor::knn(data = data, query = query, k = k,
            ...)$nn.idx
    }
    knnIdx
}

rowZscores<-function (m = NULL, min = -2, max = 2, limit = FALSE)
{
    z <- sweep(m - Matrix::rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    if (limit) {
        z[z > max] <- max
        z[z < min] <- min
    }
    return(z)
}

ArchRDORCPlot<-function(mat,
			scaleRows=TRUE,
			limits = c(-1.5, 1.5),
			rowOrder=NULL,
			labelTop=50,
			maxFeatures=nrow(mat),
			varCutOff=0.25,
			labelMarkers=NULL,
			name="DORCHeatmap",...){
    rSNA <- rowSums(is.na(mat))
    if (sum(rSNA > 0) > 0) {
        mat <- mat[rSNA == 0, ]
    }
    varQ <- getQuantiles(matrixStats::rowVars(mat))
    orderedVar <- FALSE
    if (is.null(rowOrder)) {
        mat <- mat[order(varQ, decreasing = TRUE), ]
        orderedVar <- TRUE
        if (is.null(varCutOff) & is.null(maxFeatures)) {
            n <- nrow(mat)
        }
        else if (is.null(varCutOff)) {
            n <- maxFeatures
        }
        else if (is.null(maxFeatures)) {
            n <- (1 - varCutOff) * nrow(mat)
        }
        else {
            n <- min((1 - varCutOff) * nrow(mat), maxFeatures)
        }
        n <- min(n, nrow(mat))
        mat <- mat[head(seq_len(nrow(mat)), n), ]
    }
    if (!is.null(labelTop)) {
        if (orderedVar) {
            idxLabel <- rownames(mat)[seq_len(labelTop)]
        }
        else {
            idxLabel <- rownames(mat)[order(varQ, decreasing = TRUE)][seq_len(labelTop)]
        }
    }
    else {
        idxLabel <- NULL
    }
    if (!is.null(labelMarkers)) {
        idxLabel2 <- match(tolower(labelMarkers), tolower(rownames(mat)),
            nomatch = 0)
        idxLabel2 <- idxLabel2[idxLabel2 > 0]
    }
    else {
        idxLabel2 <- NULL
    }
    idxLabel <- c(idxLabel, rownames(mat)[idxLabel2])
    if (scaleRows) {
        mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat),
            `/`)
        mat[mat > max(limits)] <- max(limits)
        mat[mat < min(limits)] <- min(limits)
    }
    if (nrow(mat) == 0) {
        stop("No Features Remaining!")
    }
    if (!is.null(rowOrder)) {
        idx <- rowOrder
    }
    else {
        idx <- order(apply(mat, 1, which.max))
    }
    #return(mat[idx,])
    p=ArchRHeatmap(mat = mat[idx, ], scale = FALSE, limits = c(min(mat),
            max(mat)), clusterCols = TRUE, clusterRows = TRUE,
            labelRows = FALSE, labelCols = TRUE, customRowLabel = match(idxLabel,
                rownames(mat[idx, ])), showColDendrogram = TRUE, draw = FALSE,name=name,...)
    return(p)
}

 
getQuantiles<-function (v = NULL, len = length(v))
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

