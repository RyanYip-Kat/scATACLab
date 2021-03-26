library(ggplot2)
library(stringr)
library(Seurat)
library(Signac)
library(argparse)

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

################################
MyplotGroups <- function(
  seurat = NULL,
  groupBy = "Sample", 
  assay = "RNA", 
  name = "Krt14",
  maxCells = 1000,
  quantCut = c(0.002, 0.998),
  pal = NULL,
  discreteSet = "stallion",
  ylim = NULL, 
  size = 0.5, 
  baseSize = 6, 
  ratioYX = NULL,
  ridgeScale = 2,
  plotAs = "violin",
  ...
  ){
  

  #Make Sure ColorBy is valid!
  if(length(assay) > 1){
    stop("assay must be of length 1!")
  }
  assayNames=names(seurat@assays)
  stopifnot(assay%in%assayNames)
  stopifnot(plotAs%in%c("ridges","violin"))

  metadata=seurat@meta.data
  groups <- metadata[,groupBy,drop=FALSE]
  groups[[groupBy]]=as.character(groups[[groupBy]])
  groupNames <- groups[,1]
  names(groupNames) <- rownames(groups)
  groupNames2 <- gtools::mixedsort(unique(groupNames))


  plotParams <- list(...)
  log2Norm <- TRUE

  Mat=GetAssayData(seurat,slot="data",assay=assay)
  colorMat=Mat[name,]

  colorList <- lapply(seq_len(nrow(colorMat)), function(x){
      colorParams <- list()
      colorParams$color <- colorMat[x, ]
      if(!is.null(discreteSet)){
        colorParams$pal <- suppressMessages(ArchR::paletteDiscrete(values = groupNames2, set = discreteSet))
      }
      if(!is.null(pal)){
        colorParams$pal <- pal
      }
      colorParams
    })


  if(!is.null(maxCells)){
    splitGroup <- split(names(groupNames), groupNames)
    useCells <- lapply(splitGroup, function(x){
      if(length(x) > maxCells){
        sample(x, maxCells)
      }else{
        x
      }
    }) %>% unlist %>% as.vector
    idx <- match(useCells, names(groupNames))
  }else{
    idx <- seq_along(groupNames)
  }

  pl <- lapply(seq_along(colorList), function(x){

    message(paste0(x, " "), appendLF = FALSE)

    if(is.null(ylim)){
      ylim <- range(colorList[[x]]$color,na.rm=TRUE) %>% extendrange(f = 0.05)
    }

    plotParamsx <- plotParams
    plotParamsx$x <- groupNames[idx]
    if(!is.null(quantCut)){
      plotParamsx$y <- quantileCut(colorList[[x]]$color[idx], min(quantCut), max(quantCut))
    }else{
      plotParamsx$y <- colorList[[x]]$color[idx]
    }
    plotParamsx$xlabel <- groupBy
    plotParamsx$ylabel <- name[x]
    plotParamsx$baseSize <- baseSize
    plotParamsx$ridgeScale <- ridgeScale
    plotParamsx$ratioYX <- ratioYX
    plotParamsx$size <- size
    plotParamsx$plotAs <- plotAs
    plotParamsx$pal <- colorList[[x]]$pal

    p <- do.call(ArchR::ggGroup, plotParamsx)

    p

  })

  names(pl) <- name
  message("")
  
  if(length(name)==1){
    pl[[1]]
  }else{
    pl
  }

}

#######################
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")


parser$add_argument("--assay",
                    type="character",
                    default="RNA",
                    choices=c("RNA","archrGA"),
                    help="gene expression assay use")

parser$add_argument("--groupby",
                    type="character",
                    default="seurat_clusters",
                    help="which column  in metadata as group")

parser$add_argument("--gene",
                    #nargs="+",
                    type="character",
                    default=NULL,
                    help="genes name to be plotted")

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
outDir=file.path(args$outdir,"VlnPlot")
makedir(outDir)

################################
message("INFO : Loading  dataset ...")
seurat=readRDS(args$seurat)
DF=read.csv(args$gene,stringsAsFactors=F,sep=",",header=F)
features=as.character(DF$V1)

message("INFO : Plot ...")
#tryCatch({
#	cat(sprintf("INFO : featurePlot [ %d ] genes\n",length(features)))
#	plotList=MyplotGroups(seurat=seurat,
#			      assay=args$assay,
#			      groupBy=args$groupby,
#			      name=features,
#			      plotAs="violin")
#	for(name in names(plotList)){
#	      	cat(sprintf("INFO : Save --- [ %s ] \n",name))
#                p=plotList[[name]]
#                ggsave(file.path(outDir,paste0(name,".pdf")),plot=p,width=12,height=10)
#	}
#
#	},error=function(e){
#		gene=features[!features%in%rownames(seurat)]
#		geneList=paste(gene,collapse=",")
#		cat(sprintf("INFO : Invalid features : [  %s  ]\n",geneList))})

cat(sprintf("INFO : featurePlot [ %d ] genes\n",length(features)))
plotList=MyplotGroups(seurat=seurat,
                              assay=args$assay,
                              groupBy=args$groupby,
                              name=features,
                              plotAs="violin")

for(name in names(plotList)){
                cat(sprintf("INFO : Save --- [ %s ] \n",name))
                p=plotList[[name]] #+ ggpubr::stat_compare_means()
                ggsave(file.path(outDir,paste0(name,".pdf")),plot=p,width=12,height=10)
        }

message("INFO : Done!")
