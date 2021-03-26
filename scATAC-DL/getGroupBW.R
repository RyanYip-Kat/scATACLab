library(ArchR)
library(Signac)
library(chromVAR)
library(motifmatchr)
library(chromVARmotifs)
library(SummarizedExperiment)
library(argparse)
parser <- ArgumentParser(description='ChromVAR handle motif from ArchR')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="ArchR Project Path")

parser$add_argument("--groupBy",
		    type="character",
		    default=NULL,
		    help="goupby use for exporting BW file in function getGroupBW")

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

outDir=args$outdir
groupBy=args$groupBy
makedir(outDir)

message("INFO : Loading dataset")
proj=loadArchRProject(args$project)
getGroupBW(proj,groupBy=groupBy)

cmd=paste0("mv ",file.path(getOutputDirectory(proj), "GroupBigWigs")," ",outDir,"/")
system(cmd)
