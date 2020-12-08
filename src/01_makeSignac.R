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
parser$add_argument("--path",
                    type="character",
                    default=NULL,
                    help="The path from cellranger-atac count or aggr")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")

args <- parser$parse_args()

###############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

make10XObject=function(path,genome){
	files=list.files(path)
	assertthat::assert_that("fragments.tsv.gz"%in%files)
	assertthat::assert_that("singlecell.csv"%in%files)
	assertthat::assert_that("filtered_peak_bc_matrix.h5"%in%files)

	countDir=file.path(path,"filtered_peak_bc_matrix.h5")
	fragDir=file.path(path,"fragments.tsv.gz")
	metaDir=file.path(path,"singlecell.csv")

	message("INFO : Create Signac Object from 10X output")
	counts <- Read10X_h5(countDir)
        metadata <- read.csv(
			     file = metaDir,
			     header = TRUE,
			     row.names = 1)
	
	seurat_assay <- CreateChromatinAssay(
					     counts = counts,
					     sep = c(":", "-"),
					     genome =genome,
					     fragments = fragDir,
					     min.features = 200,
					     min.cells = 10)
	seurat <- CreateSeuratObject(
				     counts = seurat_assay,
				     assay = 'peaks',
				     project = 'ATAC',
				     min.cells=200,
				     meta.data = metadata)

	message("INFO : Make Done!")
	return(seurat)
}


###############################  Configure
outDir=args$outdir
Genome=args$genome
Path=args$path
makedir(outDir)

if(Genome=="hg38"){
	library(EnsDb.Hsapiens.v86)
        annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
}else{
	library(EnsDb.Mmusculus.v79)
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
}

# change to UCSC style since the data was mapped to hg19
cat(sprintf("INFO : Annotation with [ %s ]\n",Genome))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- Genome


################################   Create Object  and Annotation
seurat=make10XObject(Path,Genome)
Annotation(seurat) <- annotations

############################### Computing QC Metrics
message("INFO : Computing QC Metrics")
MetricDir=file.path(outDir,"Metric")
makedir(MetricDir)


message("INFO : compute nucleosome signal score per cell")
seurat <- NucleosomeSignal(object = seurat)

message("INFO : compute TSS enrichment score per cell")
seurat <- TSSEnrichment(object = seurat, fast = FALSE)

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
############################### Save
message("INFO : Save Object")
saveRDS(seurat,file.path(outDir,"seurat.rds"))
