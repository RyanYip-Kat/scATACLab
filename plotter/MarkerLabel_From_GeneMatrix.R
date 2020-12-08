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
		    nargs="+",
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


width=args$width
height=args$height

print(paste0("# The plot size is : [ ",width,height," ]"))
markerGenes=args$markers
p <- plotEmbedding(
		   ArchRProj = projHeme,
                   colorBy = "GeneScoreMatrix",
                   name = markerGenes,
                   embedding = "UMAP",
                   quantCut = c(0.01, 0.95),
                   imputeWeights = NULL)

#p2 <- lapply(p, function(x){
#		     x + guides(color = FALSE, fill = FALSE) +
#				     theme_ArchR(baseSize = 6.5) +
#                                     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#                                     theme(
#                                          axis.text.x=element_blank(),
#                                          axis.ticks.x=element_blank(),
#                                          axis.text.y=element_blank(),
#                                          axis.ticks.y=element_blank())})

do.call(cowplot::plot_grid, c(list(ncol = 3),p))
ggsave(file.path(args$outdir,"Plot-UMAP-Marker-Genes-WO-Imputation.pdf"),width = 24,height = 24,device =cairo_pdf)
dev.off()
MyplotPDF(plotList = p,
		name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
                outpath = args$outdir,
                addDOC = FALSE, width =width, height =height)

p <- plotEmbedding(ArchRProj = projHeme,
                           colorBy = "GeneScoreMatrix",
                           name = markerGenes,
                           embedding = "UMAP",
                           quantCut = c(0.01, 0.95),
                           imputeWeights = getImputeWeights(projHeme))
do.call(cowplot::plot_grid, c(list(ncol = 3),p))
ggsave(file.path(args$outdir,"Plot-UMAP-Marker-Genes-W-Imputation.pdf"),width = 24,height = 24,device =cairo_pdf)
dev.off()
#p2 <- lapply(p, function(x){
#                             x + guides(color = FALSE, fill = FALSE) +
#                                     theme_ArchR(baseSize = 6.5) +
#                                     theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
#                                     theme(
#                                          axis.text.x=element_blank(),
#                                          axis.ticks.x=element_blank(),
#                                          axis.text.y=element_blank(),
#                                          axis.ticks.y=element_blank())
#
#                           })
MyplotPDF(plotList = p,
                name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf",
                outpath = args$outdir,
                addDOC = FALSE, width =width, height =height)

