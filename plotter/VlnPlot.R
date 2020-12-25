library(argparse)
library(stringr)
library(future)
library(Signac)
library(Seurat)
library(patchwork)
library(ggplot2)
source("/home/ye/Work/R/scATAC/Signac/plotter/plotDF.R")
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
                    nargs="+",
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
DefaultAssay(seurat)=args$assay

################################
my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

features=args$gene
plotList=list()
for(i in seq_along(features)){
	tryCatch({
		gene=features[i]
		cat(sprintf("INFO : featurePlot [ %d of %d ]  --- [ %s ]\n",i,length(features),gene))
		p=VlnPlot(object=seurat,
			      features=gene,
			      pt.size=0.1,
			      group.by=args$groupby)
			      #cols=my_palette)+theme_ArchR()
		plotList[[gene]]=p
		ggsave(file.path(outDir,paste0(gene,".pdf")),plot=p,width=12,height=10)
	},error=function(e){cat(sprintf("INFO : Invalid feature : [ %s ]\n",gene))})
}

#MyplotPDF(plotList = plotList,
#                name = "geneExpression.pdf",
#                outpath = outDir,
#                addDOC = FALSE, width =16, height = 12)



