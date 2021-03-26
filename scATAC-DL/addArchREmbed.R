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


parser$add_argument("--DF",
		    type="character",
		    default=NULL,
		    help="Embed DF csv from SCALE.py")

parser$add_argument("--name",
		    type="character",
		    default="SCALE_UMAP",
		    help="DF slot name")

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
makedir(outDir)

message("INFO : Loading dataset")
proj=loadArchRProject(args$project)
DF=read.csv(args$DF,row.names=1,stringsAsFactors=F,header=T)

message("INFO : Create Embed SimpleList")
simDF=SimpleList("df"=DF,"params"=list("SCALE"=TRUE,"n_latent"=10))
proj@embeddings[[name]]=simDF
saveRDS(proj,file.path(getOutputDirectory(proj),"Save-ArchR-Project.rds"))
