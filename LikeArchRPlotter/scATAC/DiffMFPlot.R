library(ArchR)
library(SummarizedExperiment)
library(stringr)
library(argparse)

##########################
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
  stopifnot(useName%in%names(assays(se)))
  cutOff=paste0(useName," <= 0.05 & abs(Log2FC) >= 0.5")
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
    p=ggPoint(
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
    p=ggPoint(
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
    p=ggPoint(
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
  return(p)
}

###############################
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Plot Different MF from getFeatureMarkers')
parser$add_argument("--se",
                    type="character",
                    default=NULL,
                    help="A SummarizedExperiment object returned by getMarkerFeatures()")

parser$add_argument("--useName",
                    type="character",
                    default="Pval",
                    help="useName for plotting")

parser$add_argument("--plotAs",
                    type="character",
                    default="volcano",
		    choices=c("volcano","ma","volcanodiff"),
                    help="which plot to be show")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

##############################
message("INFO : loading dataset ...")
se=readRDS(args$se)
useName=args$useName
plotAs=args$plotAs

message("INFO : Plot ...")
p=MyplotMarkers(seMarker=se,useName=useName,plotAs=plotAs)
MyplotPDF(p, name = paste0(useName,"-",plotAs,".pdf"), outpath=args$outdir, addDOC = FALSE, width = 16, height = 12)
