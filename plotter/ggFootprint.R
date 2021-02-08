library(ggplot2)
library(ArchR)

getAssay=function (se = NULL, assayName = NULL)
{
    AssayNames <- function(se) {
        names(SummarizedExperiment::assays(se))
    }
    if (is.null(assayName)) {
        o <- SummarizedExperiment::assay(se)
    }
    else if (assayName %in% AssayNames(se)) {
        o <- SummarizedExperiment::assays(se)[[assayName]]
    }
    else {
        stop(sprintf("assayName '%s' is not in assayNames of se : %s",
            assayName, paste(AssayNames(se), collapse = ", ")))
    }
    return(o)
}

ggFootprint <- function(
  seFoot = NULL,
  name = NULL,
  pal = NULL,
  smoothWindow = 6,
  flank = 250,
  flankNorm = 50,
  baseSize = 6,
  normMethod = "subtract",
  cutoffX=100
  ){

  errorList <- list()

  #Get Footprint Info
  rowDF <- SummarizedExperiment::rowData(seFoot)
  footMat <- getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="footprint"),], name)
  biasMat <- getAssay(seFoot[BiocGenerics::which(rowDF[,2]=="bias"),], name)
  footDF <- rowDF[BiocGenerics::which(rowDF[,2]=="footprint"),]
  biasDF <- rowDF[BiocGenerics::which(rowDF[,2]=="bias"),]

  errorList$footMat <- footMat
  errorList$biasMat <- biasMat
  errorList$footDF <- footDF
  errorList$biasDF <- biasDF

  #Smooth Foot and Bias Mat because of sparsity
  if(!is.null(smoothWindow)){
    message("Applying smoothing window to footprint ...")
    footMat <- apply(footMat, 2, function(x) ArchR:::.centerRollMean(x, smoothWindow))
    biasMat <- apply(biasMat, 2, function(x) ArchR:::.centerRollMean(x, smoothWindow))
  }

  #Normalize Foot and Bias Mat
  DFF=flank - flankNorm
  idx <- which(abs(footDF$x) >= DFF)
  footMat <- t(t(footMat) / colMeans(footMat[idx, ,drop=FALSE]))
  biasMat <- t(t(biasMat) / colMeans(biasMat[idx, ,drop=FALSE]))

  errorList$footMatNorm <- footMat
  errorList$biasMatNorm <- footMat

  #Norm Foot By Bias
  if(tolower(normMethod) == "none"){
    title <- ""
  }else if(tolower(normMethod) == "subtract"){
    title <- "Tn5 Bias Subtracted\n"
    footMat <- footMat - biasMat
  }else if(tolower(normMethod) == "divide"){
    title <- "Tn5 Bias Divided\n"
    footMat <- footMat / biasMat
  }else{
    stop("normMethod not recognized!")
  }

  #Get Mean and SD for each Assay
  footMatMean <- ArchR:::.groupMeans(footMat, SummarizedExperiment::colData(seFoot)$Group)
  footMatSd <- ArchR:::.groupSds(footMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatMean <- ArchR:::.groupMeans(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  biasMatSd <- ArchR:::.groupSds(biasMat, SummarizedExperiment::colData(seFoot)$Group)
  smoothFoot <- rowMaxs(apply(footMatMean, 2, function(x) ArchR:::.centerRollMean(x, 11)))

  errorList$footMatMean <- footMatMean
  errorList$footMatSd <- footMatSd
  errorList$biasMatMean <- biasMatMean
  errorList$biasMatSd <- biasMatSd
  errorList$smoothFoot <- smoothFoot

  #Create Plot Data Frames
  plotIdx <- seq_len(nrow(footMatMean)) #sort(unique(c(1, seq(1, nrow(footMatMean), smoothWindow), nrow(footMatMean))))
  plotFootDF <- lapply(seq_len(ncol(footMatMean)), function(x){
    data.frame(
      x = footDF$x, 
      mean = footMatMean[,x], 
      sd = footMatSd[,x], 
      group = colnames(footMatMean)[x]
      )[plotIdx,,drop=FALSE]
  }) %>% Reduce("rbind",. )
  plotFootDF$group <- factor(paste0(plotFootDF$group), levels = unique(gtools::mixedsort(paste0(plotFootDF$group))))

  plotBiasDF <- lapply(seq_len(ncol(biasMatMean)), function(x){
    data.frame(
      x = biasDF$x, 
      mean = biasMatMean[,x], 
      sd = biasMatSd[,x], 
      group = colnames(biasMatMean)[x]
      )[plotIdx,,drop=FALSE]
  }) %>% Reduce("rbind",. )
  plotBiasDF$group <- factor(paste0(plotBiasDF$group), levels = unique(gtools::mixedsort(paste0(plotBiasDF$group))))

  errorList$plotFootDF <- plotFootDF
  errorList$plotBiasDF <- plotBiasDF

  out <- tryCatch({

    #Plot GG
    if(is.null(pal)){
      pal <- paletteDiscrete(values=gtools::mixedsort(SummarizedExperiment::colData(seFoot)$Group))
    }

    
    if(!is.null(cutoffX)){
	    plotFootDF=subset(plotFootDF,abs(x) < cutoffX)
    }
    plotFootDF$mean=  100*exp(plotFootDF$mean) #-log2(plotFootDF$mean+max(abs(plotFootDF$mean)))
    plotFootDF$sd= 10*exp(plotFootDF$sd) #-log2(plotFootDF$sd+max(abs(plotFootDF$sd)))

    plotMax <- plotFootDF[order(plotFootDF$mean,decreasing=TRUE),]
    plotMax <- plotMax[abs(plotMax$x) > 20 & abs(plotMax$x) < 50, ] #<= flank - flankNorm,]
    plotMax <- plotMax[!duplicated(plotMax$group),]
    plotMax <- plotMax[seq_len(ceiling(nrow(plotMax) / 4)), ]
    plotMax$x <- 25


    ggFoot <- ggplot(plotFootDF, aes(x = x, y = mean, color = group)) + 
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
      geom_line() + 
      scale_color_manual(values = pal) + 
      scale_fill_manual(values = pal) + 
      xlab("Distance to motif center (bp)") +
      coord_cartesian(
        expand = FALSE, 
        #ylim = c(quantile(plotFootDF$mean, 0.0001), 1.15*quantile(smoothFoot, 0.999)), 
        xlim = c(min(plotFootDF$x),max(plotFootDF$x))
      ) + theme_ArchR(baseSize = baseSize) + ggtitle(name) +
      guides(fill = FALSE) + 
      guides(color = FALSE) + ylab(paste0(title,"Normalized Insertions")) +
      ggrepel::geom_label_repel(data = plotMax, aes(label = group), size = 3, xlim = c(75, NA))

    if(!is.null(cutoffX)){
            plotBiasDF=subset(plotBiasDF,abs(x) < cutoffX)
    }
    plotBiasDF$mean=100*exp(plotBiasDF$mean) #-log2(plotBiasDF$mean+max(abs(plotBiasDF$mean)))
    plotBiasDF$sd=10*exp(plotBiasDF$sd) #-log2(plotBiasDF$sd+max(abs(plotBiasDF$sd)))

    ggBias <- ggplot(plotBiasDF, aes(x = x, y = mean, color = group)) + 
      geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd, linetype = NA, fill = group), alpha = 0.4) +
      geom_line() + 
      scale_color_manual(values = pal) + 
      scale_fill_manual(values = pal) + 
      xlab("Distance to motif center (bp)") +
      coord_cartesian(
        expand = FALSE, 
        #ylim = c(quantile(plotBiasDF$mean, 0.0001), 1.05*quantile(plotBiasDF$mean, 0.999)), 
        xlim = c(min(plotBiasDF$x),max(plotBiasDF$x))
      ) + theme_ArchR(baseSize = baseSize) + ylab("Tn5-Bias Normalized Insertions") + 
      theme(legend.position = "bottom", legend.box.background = element_rect(color = NA)) 
    list("FootDF"=plotFootDF,"BiasDF"=plotBiasDF,"ggFoot"=ggFoot,"ggBias"=ggBias)       
    #ggAlignPlots(ggFoot, ggSmallLegend(ggBias), sizes=c(2,1), draw = FALSE)

  }, error = function(e){message("INFO : Plot failed!")
  })

  out

}

ggSmallLegend <- function(
  gg = NULL,
  pointSize = 2,
  baseSize = 5,
  spaceLegend = 0.1
  ) {
    #https://stackoverflow.com/questions/52297978/decrease-overal-legend-size-elements-and-text
    gg +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
               color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = baseSize), 
              legend.text  = element_text(size = baseSize),
              legend.key.size = unit(spaceLegend, "lines"))
}
