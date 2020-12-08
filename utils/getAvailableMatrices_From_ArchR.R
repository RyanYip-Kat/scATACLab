library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

parser$add_argument("--slot",
		    type="character",
		    default="GeneScoreMatrix")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print(paste0("### Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("### Loading Project")
projHeme<-loadArchRProject(args$project)

avaliable_matrixs=getAvailableMatrices(projHeme)
stopifnot(args$slot%in%avaliable_matrixs)
se=getMatrixFromProject(projHeme,useMatrix=args$slot)

feature=getPeakSet(projHeme)
if(args$slot=="GeneScoreMatrix"){
	features=getFeatures(projHeme,useMatrix =args$slot)
        rownames(se)=features
}else if(args$slot=="PeakMatrix"){
	peaks=getPeakSet(projHeme)
	# mcols(peaks)
	rowData(se)=cbind(rowData(se),mcols(peaks))
	rgs=as.data.frame(ranges(peaks))
	sqs=as.data.frame(seqnames(peaks))
	peakset=cbind(sqs,rgs)
	peak_set=paste(peakset$value,peakset$start,peakset$end,sep="_")
	rownames(se)=peak_set
}

print(paste0("### Save ",args$slot," Matrix"))
saveRDS(se,file.path(args$outdir,paste0(args$slot,"_Summarized-Experiment.rds")))
