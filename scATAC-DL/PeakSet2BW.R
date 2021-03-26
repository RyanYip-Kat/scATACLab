library(ArchR)
library(rtracklayer)
library(argparse)
parser <- ArgumentParser(description='ChromVAR handle motif from ArchR')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="ArchR Project Path")


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
peakSet=getPeakSet(ArchRProj=proj)
len=getChromLengths(proj)[names(seqlengths(peakSet))]
seqlengths(peakSet)=len

message("INFO : Export Into BW File")
filename=file.path(outDir,"PeakSet.bw")
rtracklayer::export(peakSet,con=filename,format="BigWig")
message("INFO : Done!")
