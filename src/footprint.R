library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="the path of project saved")

parser$add_argument("--motif",
                    type="character",
                    default=NULL,
                    help="motif txt or csv file:only one column")

parser$add_argument("--group",
                    type="character",
                    default="Clusters",
                    help="the group by  name")

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

width=12
height=12

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads)

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)

print("# Get Motif position")
motifPositions <- getPositions(projHeme)

motifs=read.table(args$motif,sep="\t",stringsAsFactors=FALSE,header=FALSE)$V1
#motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]

meta=as.data.frame(getCellColData(projHeme))
stopifnot(args$group%in%colnames(meta))
print("# getFootprints")
seFoot <- getFootprints(
  ArchRProj = projHeme,
  positions = motifPositions[markerMotifs],
  groupBy =args$group
)


saveArchRProject(ArchRProj = projHeme, outputDirectory = file.path(args$outdir,"Save-ProjHeme-footprint"), load = TRUE)
