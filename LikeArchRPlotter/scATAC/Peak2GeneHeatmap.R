# ArchR > 1.0.0
library(argparse)
library(stringr)
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


parser$add_argument("--groupby",
                    type="character",
                    default="Clusters",
                    help="the project path  of ArchR")

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

width=12
height=10

message("INFO : Loading dataset ...")
projHeme=loadArchRProject(args$project)
peakSet=getPeakSet(projHeme)

message("INFO : Plot P2G heatmap ...")
p=plotPeak2GeneHeatmap(projHeme,
		       groupBy=args$groupby,
		       corCutOff = 0.45,
		       FDRCutOff = 1e-04,
		       k = 15,
		       nPlot = 25000,
		       limitsATAC = c(-2, 2),
		       limitsRNA = c(-2, 2),
		       palGroup = NULL,
		       palATAC = paletteContinuous("solarExtra"),
		       palRNA = paletteContinuous("blueYellow"),
		       verbose = TRUE,
		       returnMatrices=FALSE)

MyplotPDF(p, name = paste0("P2GHeatmap-",args$groupby,".pdf"), outpath=args$outdir, addDOC = FALSE, width = width, height = height)
pmat=plotPeak2GeneHeatmap(projHeme,
                       groupBy=args$groupby,
                       corCutOff = 0.45,
                       FDRCutOff = 1e-04,
                       k = 15,
                       nPlot = 25000,
                       limitsATAC = c(-2, 2),
                       limitsRNA = c(-2, 2),
                       palGroup = NULL,
                       palATAC = paletteContinuous("solarExtra"),
                       palRNA = paletteContinuous("blueYellow"),
                       verbose = TRUE,
                       returnMatrices=TRUE)

############################
Peak2GeneLinks=pmat[["Peak2GeneLinks"]]
ATAC=pmat[["ATAC"]]
RNA=pmat[["RNA"]]

kmeansId=ATAC[["kmeansId"]]
DF=data.frame("peak"=Peak2GeneLinks$peak,"gene"=Peak2GeneLinks$gene,"kid"=as.factor(kmeansId))
rownames(DF)=rownames(Peak2GeneLinks)

N=100
DF_list=list()
for(id in unique(DF$kid)){
        cat(sprintf("INFO : select [ %d ] gene in [ %s ] kmeansId\n",N,id))
        df=subset(DF,kid==id)
        n=N
        if(N>nrow(df)){
                   n=nrow(df)
        }
        x=df[1:n,,drop=FALSE]
        x$rowName=rownames(x)
        DF_list[[id]]=x
}

topDF=do.call(rbind,DF_list)
rownames(topDF)=topDF$rowName

topDF$chr=unlist(lapply(topDF$peak,function(x)str_split(x,":")[[1]][1]))
topDF$start=unlist(lapply(topDF$peak,function(x)str_split(x,":|-")[[1]][2]))
topDF$end=unlist(lapply(topDF$peak,function(x)str_split(x,":|-")[[1]][3]))
topDFGR=makeGRangesFromDataFrame(topDF,keep.extra.columns=T)


###########################
idx=findOverlaps(topDFGR,peakSet)
topDFGR=topDFGR[idx@from,]
mcols(topDFGR)$nearstGene=mcols(peakSet[idx@to,])$nearestGene
saveRDS(topDFGR,file.path(args$outdir,"topDFGR.rds"))
saveRDS(pmat,file.path(args$outdir,"Peak2GeneHeatmapMatrix.rds"))

df=mcols(topDFGR)[,c("gene","peak","nearstGene","kid")]
write.table(df,file.path(args$outdir,paste0("top_",N,".csv")),quote=F,row.names=F,sep=",")
message("INFO : Done!")

