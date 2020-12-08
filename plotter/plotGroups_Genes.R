library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
library(Cairo)
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the project path  of ArchR")


parser$add_argument("--markers",
		    #nargs="+",
                    type="character",
                    default=NULL)

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--group",
                    type="character",
                    default="Clusters")

parser$add_argument("--width",
                    type="integer",
                    default=12,
                    help="the width of plot")

parser$add_argument("--height",
                    type="integer",
                    default=12,
                    help="the height of plot")

parser$add_argument("-slot",
                    type="character",
                    default="GeneScoreMatrix",
                    help="GeneScoreMatrix,MotifMatrix")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
projHeme=loadArchRProject(args$project)


width=args$width
height=args$height

print(paste0("# The plot size is : [ ",width,height," ]"))
markers=as.character(read.table(args$markers,stringsAsFactors=F,header=F)$V1)
for(gene in markers){
	print(paste0("--- Plot :",gene))
	plotGroups(ArchRProj=projHeme,
		     colorBy=args$slot,
		     groupBy=args$group,
		     name=gene,
		     plotAs = "ridges")
	filename=file.path(args$outdir,paste0(gene,"_ridges.pdf"))
	ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
        dev.off()


	plotGroups(ArchRProj=projHeme,
                     colorBy=args$slot,
                     groupBy=args$group,
		     name=gene,
		     groupOrder=NULL,
                     plotAs = "violin")
        filename=file.path(args$outdir,paste0(gene,"_violin.pdf"))
        ggsave(filename = filename,device = cairo_pdf,height = height,width=width)
        dev.off()
}

