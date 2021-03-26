library(argparse)
library(stringr)
library(future)
library(Signac)
library(Seurat)
library(patchwork)
library(ggplot2)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")



parser$add_argument("--motif",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="motif name to be plotted")

parser$add_argument("--groupby",
                    type="character",
                    default="seurat_clusters",
                    help="group use for CoveragePlot")

parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")

parser$add_argument("--assay",
                    type="character",
                    default="peaks",
		    choices=c("peaks","ATAC","atac","chromvar"),
		    help="assay to use")

parser$add_argument("--outdir",
                    type="character",
                    default="./Results")

args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

################################
outDir=file.path(args$outdir,"footprintPlot")
genome=args$genome
makedir(outDir)

##############################
if(genome=="hg38"){
        library(BSgenome.Hsapiens.UCSC.hg38)
        BS=BSgenome.Hsapiens.UCSC.hg38
        species="Homo sapiens"
}else{
        library(BSgenome.Mmusculus.UCSC.mm10)
        BS=BSgenome.Mmusculus.UCSC.mm10
        species="Mus musculus"
}

################################
message("INFO : Loading  dataset ...")
seurat=readRDS(args$seurat)
assayNames=Assays(seurat)
#assertthat::assert_that("motifs"%in%assayNames)

################################ motifs footprint plot
motif.names=args$motif
cat(sprintf("INFO : Input [ %d ] motifs to be plotted\n",length(motif.names)))

message("INFO : gather the footprinting information for sets of motifs")
seurat <- Footprint(
  object = seurat,
  motif.name =motif.names,
  genome = BS,
  assay=args$assay
)

for(motif in motif.names){
	cat(sprintf("INFO : Plot motif : [ %s of %d]\n",motif,length(motif.names)))
	tryCatch({ p=PlotFootprint(object = seurat,
                       features=motif,
		       label.top=3,
		       assay=args$assay,
                       group.by=args$groupby)+patchwork::plot_layout(ncol = 1)
                   fileName=file.path(outDir,paste0(motif,"-motifs.pdf"))
                   ggsave(fileName,plot=p,width=8,height=16)},
                   error=function(e){cat(sprintf("INFO : Invalid Motif\n"))},
                   finnally={cat(sprintf("INFO : footprint plot [ %s ] sucessfuly!\n",motif))}
        )
}

