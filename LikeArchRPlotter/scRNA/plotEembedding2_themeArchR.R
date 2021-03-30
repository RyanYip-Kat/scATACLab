library(Seurat)
library(Signac)
library(ggplot2)
library(stringr)

source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
source("/home/ye/Work/R/scATAC/Signac/plotter/SeuratImputeWeight.R")

tolower<-str_to_lower
quantileCut<-function (x = NULL, lo = 0.025, hi = 0.975, maxIf0 = TRUE)
{
    q <- quantile(x, probs = c(lo, hi),na.rm=TRUE)
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

MyplotEmbedding <- function(
  seurat = NULL,
  assay="RNA",
  embedding = "UMAP",
  colorBy = "metadata", # column in metadata or genes in matrixs
  name = "Sample",
  imputeWeights=NULL,
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
  embedding=tolower(embedding)
  stopifnot(embedding%in%Reductions(seurat))
  stopifnot(colorBy%in%c("metadata","matrix"))  # only support two style
  df <- Embeddings(seurat,reduction=embedding)
  meta.data=seurat@meta.data
  if(!all(rownames(df) %in% colnames(seurat))){
    stop("Not all cells in embedding are present in ArchRProject!")
  }

  if(!is.null(sampleCells)){
    if(sampleCells < nrow(df)){
      df <- df[sort(sample(seq_len(nrow(df)), sampleCells)), , drop = FALSE]
    }
  }

  #Parameters
  plotParams <- list(...)
  plotParams$x <- df[,1]
  plotParams$y <- df[,2]
  plotParams$title <- paste0(embedding, " of ", Project(seurat))
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

  #Make Sure ColorBy is valid!
  if(length(colorBy) > 1){
    stop("colorBy must be of length 1!")
  }

  if(tolower(colorBy) == "metadata"){
      
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
  }else if(tolower(colorBy) == "matrix"){
    Mat=GetAssayData(seurat,slot="data",assay=assay)
    colorMat=Mat[name,,drop=FALSE]

    if(!all(rownames(df) %in% colnames(colorMat))){
      stop("Not all cells in embedding are present in matrix. This may be due to using a custom embedding.")
    }

    colorMat <- colorMat[,rownames(df), drop=FALSE]
    if(!is.null(imputeWeights)){
      message("Imputing Matrix")
      colorMat <- MyimputeMatrix(mat = colorMat, imputeWeights = imputeWeights)
      if(!inherits(colorMat, "matrix")){
        colorMat <- matrix(colorMat, ncol = nrow(df))
        colnames(colorMat) <- rownames(df)
      }
    }
    colorList <- lapply(seq_len(nrow(colorMat)), function(x){
      colorParams <- list()
      colorParams$color <- colorMat[x, ]
      colorParams$discrete <- FALSE
      colorParams$title <- sprintf("%s colored by\n%s : %s", plotParams$title, colorBy, name[x])
      if(tolower(colorBy) == "matrix"){
        colorParams$continuousSet <- "horizonExtra"
      }else{
        colorParams$continuousSet <- "solarExtra"
      }
      if(!is.null(continuousSet)){
        colorParams$continuousSet <- continuousSet
      }
      if(!is.null(discreteSet)){
        colorParams$discreteSet <- discreteSet
      }
      return(colorParams)
    })

  }else{
	  stop("Invalid colorBy!!!")
  }
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
MyPlot<-function(seurat,
		 assay,
		 colorBy,
		 features,
		 embedding,
		 outDir,
		 imputeWeights=NULL,
		 splitBy=NULL,
		 combine=FALSE){
	if(!is.null(splitBy)){
		seurat_list=SplitObject(seurat,split.by=splitBy)
	}else{
		seurat_list=list("SeuratObject"=seurat)  # only one object list
	}

	for(i in seq_along(features)){
		feature=features[i]
		cat(sprintf("INFO : [ %d of %d ] --- [ %s ]\n",i,length(features),feature))
		plotList=lapply(names(seurat_list),function(i){
				sr=seurat_list[[i]]
				p=MyplotEmbedding(seurat=sr,
						  assay=assay,
						  embedding=embedding,
                                                  name=feature,
						  imputeWeights=imputeWeights,
                                                  colorBy=colorBy)+ggtitle(paste0(i,"-",feature))
				rm(sr)
				return(p)
		 })
		names(plotList)=names(seurat_list)
		if(combine){
			MyplotPDF(plotList,name=feature,outpath=outDir,width=16,height=16)
		}else{
			for(name in names(plotList)){
				cat(sprintf("INFO : Save --- [ %s ] \n",name))
                                p=plotList[[name]] #+ ggpubr::stat_compare_means()
                                if(args$colorBy=="matrix"){
					ggsave(file.path(outDir,paste0(name,"-",feature,".pdf")),plot=p,width=12,height=10)
				}else{
                                        svg(file.path(outDir,paste0(name,"-",feature,".svg")),width=12,height=10)
                                        print(p)
                                        dev.off()
				}
			}
		}
	}
}


###############################
library(argparse)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
                    help="the dataset txt file,2 columns with \t")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="")

parser$add_argument("--assay",
                    type="character",
                    default="RNA",
                    choices=c("RNA","archrGA"),
                    help="gene expression assay use")

parser$add_argument("--colorBy",
                    type="character",
                    default="metadata",
		    choices=c("metadata","matrix"),
                    help="metadata or matrix")

parser$add_argument("--splitBy",
		    type="character",
		    default=NULL,
		    help="split by column in metadata")

parser$add_argument("--name",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="genes name or column to be plotted")

parser$add_argument("--embed",
                    type="character",
                    default="umap")

parser$add_argument("--outdir",
                    type="character",
                    default="output")
args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

################################
outDir=file.path(args$outdir,"Embedding")
makedir(outDir)

################################
message("INFO : Loading  dataset ...")
seurat=readRDS(args$seurat)
if(!is.null(args$column) & !is.null(args$subset)){
        #Idents(seurat)<-seurat$status
        #seurat<-subset(seurat,idents=args$status)
        metadata=seurat@meta.data
        Idents(seurat)=metadata[[args$column]]
        seurat<-subset(seurat,idents=args$subset)
}

if(args$colorBy=="matrix"){
	message("INFO : Impute Weights")
        impWeights=ImputeWeights(seurat,
			 reducedDims="pca",
			 sampleCells=10000,
			 nRep = 2,
			 threads=8)
}else{
	impWeights=NULL
}

#DF=read.csv(args$name,stringsAsFactors=F,sep=",",header=F)
#features=as.character(DF$V1)
features=args$name
message("INFO : Plot ...")
MyPlot(seurat=seurat,
       assay=args$assay,
       splitBy=args$splitBy,
       embedding=args$embed,
       imputeWeights=impWeights,
       features=features,
       colorBy=args$colorBy,
       outDir=outDir,
       combine=FALSE)
message("INFO : Done!")

