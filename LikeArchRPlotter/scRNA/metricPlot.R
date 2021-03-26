library(stringr)
library(DropletUtils)
library(ggplot2)
library(Seurat)
library(Signac)
library(argparse)
library(GenomicRanges)
library(future)
#############################
source("/home/ye/Work/R/scATAC/Signac/src/helper_functions.R")
parser <- ArgumentParser(description='Program to make Signac Object')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="seurat object from make Signac")


args <- parser$parse_args()

###############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

seurat=readRDS(args$seurat)
MetricDir=file.path(dirname(args$seurat),"Metric")
makedir(MetricDir)


message("INFO : add blacklist ratio and fraction of reads in peaks")
seurat$pct_reads_in_peaks <- seurat$peak_region_fragments / seurat$passed_filters * 100
seurat$blacklist_ratio <- seurat$blacklist_region_fragments / seurat$peak_region_fragments

seurat$high.tss <- ifelse(seurat$TSS.enrichment > 2, 'High', 'Low')
p=TSSPlot(seurat, group.by = 'high.tss') + NoLegend()
ggsave(file.path(MetricDir,"high.tss.pdf"),plot=p,width=10,height=12)

seurat$nucleosome_group <- ifelse(seurat$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p=FragmentHistogram(object = seurat, group.by = 'nucleosome_group')
ggsave(file.path(MetricDir,"nucleosome_group.pdf"),plot=p,width=10,height=12)

p=VlnPlot(
  object = seurat,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

ggsave(file.path(MetricDir,"multiple_metric.pdf"),plot=p,width=24,height=10)
message("INFO : Sample Metric Done!")

