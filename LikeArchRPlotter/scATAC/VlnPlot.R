library(argparse)
library(stringr)
library(ggplot2)
library(ArchR)

parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="ArchR Project Path")


parser$add_argument("--groupby",
                    type="character",
                    default="seurat_clusters",
                    help="which column  in metadata as group")

parser$add_argument("--gene",
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
addArchRThreads(threads =8)

message("INFO : Loading dataset ...")
projHeme=loadArchRProject(args$project)
DF=read.csv(args$gene,stringsAsFactors=F,sep=",",header=F)
features=as.character(DF$V1)

plotList=plotGroups(ArchRProj=projHeme,
			     groupBy=args$groupby,
			     colorBy="GeneScoreMatrix",
			     name=features,
			     plotAs="violin")

for(i in seq_along(plotList)){
        name=names(plotList)[i]
        cat(sprintf("INFO : Save --- [ %d of %d ] --- [ %s ]\n",i,length(plotList),name))
        p=plotList[[i]]
        ggsave(file.path(args$outdir,paste0(name,".pdf")),plot=p,width=12,height=10)
}

message("INFO : Done")

