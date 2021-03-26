library(ggplot2)
library(stringr)

source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
tolower<-str_to_lower
quantileCut<-function (x = NULL, lo = 0.025, hi = 0.975, maxIf0 = TRUE)
{
    q <- quantile(x, probs = c(lo, hi))
    if (q[2] == 0) {
        if (maxIf0) {
            q[2] <- max(x)
        }
    }
    x[x < q[1]] <- q[1]
    x[x > q[2]] <- q[2]
    return(x)
}

isDiscrete<-function (x = NULL)
{
    is.factor(x) || is.character(x) || is.logical(x)
}


mergeParams<-function (paramInput = NULL, paramDefault = NULL)
{
    for (i in seq_along(paramDefault)) {
        if (!(names(paramDefault)[i] %in% names(paramInput))) {
            paramInput[[names(paramDefault)[i]]] <- paramDefault[[i]]
        }
    }
    return(paramInput)
}

plotSCALE <- function(
  emb, # UMAP or tSNE from SCALE
  metadata,  # meta data for each Cell
  name = "Clusters",  # colorBy by name
  pal = NULL,
  size = 0.1,
  sampleCells = NULL,
  highlightCells = NULL,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10,
  plotAs = NULL,
  labelSize=0, #no label
  ...
  ){


  ##############################
  # Get Embedding
  ##############################
  #df <- getEmbedding(ArchRProj, embedding = embedding, returnDF = TRUE)
  df <- emb
  meta.data=metadata
  meta.data=meta.data[rownames(df),]

  if(!is.null(sampleCells)){
    if(sampleCells < nrow(df)){
      df <- df[sort(sample(seq_len(nrow(df)), sampleCells)), , drop = FALSE]
    }
  }

  #Parameters
  plotParams <- list(...)
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- "SCALE of scATAC"
  plotParams$baseSize <- baseSize
  
  #Additional Params!
  plotParams$xlabel <- colnames(df)[1]
  plotParams$ylabel <- colnames(df)[2]
  plotParams$rastr <- rastr
  plotParams$size <- size
  plotParams$randomize <- randomize

  #Check if Cells To Be Highlighed
  if(!is.null(highlightCells)){
    highlightPoints <- match(highlightCells, rownames(df), nomatch = 0)
    if(any(highlightPoints==0)){
      stop("highlightCells contain cells not in Embedding cellNames! Please make sure that these match!")
    }
  }

      
  colorList <- lapply(seq_along(name), function(x){
      colorParams <- list()
      #colorParams$color <- as.vector(getCellColData(ArchRProj, select = name[x], drop = FALSE)[rownames(df), 1])
      colorParams$color <- as.vector(subset(meta.data,select=name[x])[rownames(df), 1])
      colorParams$discrete <- isDiscrete(colorParams$color)
      colorParams$continuousSet <- "solarExtra"
      colorParams$discreteSet <- "stallion"
      colorParams$title <- paste(plotParams$title, " colored by\ncolData : ", name[x])
      if(!is.null(continuousSet)){
        colorParams$continuousSet <- continuousSet
      }
      if(!is.null(discreteSet)){
        colorParams$discreteSet <- discreteSet
      }
      return(colorParams)
  })
  message("Plotting Embedding")

  ggList <- lapply(seq_along(colorList), function(x){

    message(x, " ", appendLF = FALSE)

    plotParamsx <- mergeParams(colorList[[x]], plotParams)
    if(plotParamsx$discrete){
      plotParamsx$color <- paste0(plotParamsx$color)
    }

    if(!plotParamsx$discrete){

      plotParamsx$color <- quantileCut(plotParamsx$color, min(quantCut), max(quantCut))

      plotParamsx$pal <- ArchR::paletteContinuous(set = plotParamsx$continuousSet)

      if(!is.null(pal)){

        plotParamsx$pal <- pal
        
      }

      if(is.null(plotAs)){
        plotAs <- "hexplot"
      }

      if(tolower(plotAs) == "hex" | tolower(plotAs) == "hexplot"){

        plotParamsx$discrete <- NULL
        plotParamsx$continuousSet <- NULL
        plotParamsx$rastr <- NULL
        plotParamsx$size <- NULL
        plotParamsx$randomize <- NULL

        gg <- do.call(ArchR::ggHex, plotParamsx)

      }else{

        if(!is.null(highlightCells)){
          plotParamsx$highlightPoints <- highlightPoints
        }

        gg <- do.call(ArchR::ggPoint, plotParamsx)

      }

    }else{
      
      if(!is.null(pal)){
        plotParamsx$pal <- pal
      }

      if(!is.null(highlightCells)){
        plotParamsx$highlightPoints <- highlightPoints
      }

      gg <- do.call(ArchR::ggPoint, plotParamsx)

    }

    if(!keepAxis){
      gg <- gg + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }

    gg

  })
  names(ggList) <- name
  message("")

  if(length(ggList) == 1){
    ggList <- ggList[[1]]
  }
  ggList

}

##############################
library(argparse)
parser <- ArgumentParser(description='Plot Embedding')
parser$add_argument("--embed",
                    type="character",
                    default=NULL,
                    help="Emebdding csv file")


parser$add_argument("--metadata",
                    type="character",
                    default=NULL,
                    help="meta data csv file")

parser$add_argument("--groupBy",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="group By plot list")


parser$add_argument("--outdir",
                    type="character",
                    default="output")
args <- parser$parse_args()

makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

outDir=args$outdir
makedir(outDir)
###########################
message("INFO : Loading  dataset ")
DF=read.csv(args$embed,row.names=1,stringsAsFactors=F,header=T)
metadata=read.csv(args$metadata,stringsAsFactors=F,header=T)
groups=args$groupBy
plotLists=plotSCALE(emb=DF,metadata=metadata,name=groups)
MyplotPDF(plotLists,name="SCALE-Embed",width=12,height=10,outpath=outDir)


