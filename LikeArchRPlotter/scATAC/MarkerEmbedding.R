library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the project path  of ArchR")


parser$add_argument("--marker",
		    #nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--width",
                    type="integer",
                    default=12,
                    help="the width of plot")

parser$add_argument("--height",
                    type="integer",
                    default=12,
                    help="the height of plot")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
projHeme=loadArchRProject(args$project)
if(is.null(getImputeWeights(projHeme))){
	projHeme=addImputeWeights(projHeme)
}
width=args$width
height=args$height

DF=read.csv(args$marker,stringsAsFactors=F,sep=",",header=F)
features=as.character(DF$V1)
all.genes=getFeatures(projHeme,useMatrix="GeneScoreMatrix")
features=features[features%in%all.genes]
pL <- plotEmbedding(ArchRProj = projHeme,
                   colorBy = "GeneScoreMatrix",
                   name = features,
                   embedding = "UMAP",
                   quantCut = c(0.01, 0.95),
                   imputeWeights = getImputeWeights(projHeme))

MyplotPDF(plotList =pL,
                name = "GeneScoreMatrix-Impute",
                outpath = args$outdir,
                addDOC = FALSE, width =width, height =height)


for(i in seq_along(pL)){
	name=names(pL)[i]
	cat(sprintf("INFO : Save --- [ %d of %d ] --- [ %s ]\n",i,length(pL),name))
	p=pL[[i]]
	ggsave(file.path(args$outdir,paste0(name,".pdf")),plot=p,width=width,height=height)
}



