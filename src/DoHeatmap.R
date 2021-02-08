library(argparse)
library(stringr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(future)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
		    default=NULL,
                    help="the path of project saved")

parser$add_argument("--outdir",
                    type="character",
                    default="output")

parser$add_argument("--groupby",
                    type="character",
                    default="Clusters")

parser$add_argument("--jobs",
                    type="integer",
                    default=16,
                    help="number of threads for command")

parser$add_argument("--n_genes",
                    type="integer",
                    default=10,
                    help="number of top genes to plot in heatmap")

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

#plan("multiprocess", workers = args$jobs)
print("### Loading data")
seurat=readRDS(args$seurat)

if(nrow(seurat[["RNA"]]@scale.data)==0){
	print("### Scale data")
	seurat <- ScaleData(seurat, features=VariableFeatures(object =seurat))
}
metadata=seurat@meta.data
stopifnot(args$groupby%in%colnames(metadata))
Idents(seurat)=metadata[[args$groupby]]

print("### Find Markers")
markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,file.path(args$outdir,"FindAllMarkers.csv"),sep=",",quote=FALSE,row.names=FALSE)

top_markers <- markers %>% group_by(cluster) %>% top_n(n = args$n_genes, wt = avg_logFC)
write.table(top_markers,file.path(args$outdir,paste0("FindAllMarkers_TOP_",args$n_genes,".csv")),sep=",",quote=FALSE,row.names=FALSE)


DoHeatmap(seurat, features = top_markers$gene) + NoLegend()
ggsave(file.path(args$outdir,paste0("DoHeatmap_TOP_",args$n_genes,".pdf")),width = 24,height = 24,device =cairo_pdf)
dev.off()
