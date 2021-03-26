library(argparse)
library(stringr)
library(ArchR)
parser <- ArgumentParser(description='Export BED file from reproduciblePeaks.gr.rds after call peaks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="the project path  of ArchR,must include PeakCalls subdir")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

stopifnot("PeakCalls"%in%list.files(args$project))
path=file.path(args$project,"PeakCalls")
reproduciblePeaks=list.files(path,"*reproduciblePeaks.gr.rds")

for(file in reproduciblePeaks){
	group=str_split(file,"-")[[1]][1]
	group=str_replace(group,"\\.","-")
	file=file.path(path,file)
	cat(sprintf("### Loading Peaks GR file from : %s\n",file))
	gr=readRDS(file)

	out=as.data.frame(gr)
        out1=out[,c("seqnames","start","end","GroupReplicate","score")]

	outname=file.path(args$outdir,paste0("macs2-",group,".bed"))
	cat(sprintf("### Write Peaks BED file into : %s\n",outname))
	write.table(out1,outname,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)

	outname=file.path(args$outdir,paste0("macs2-",group,"-gr.bed"))
	write.table(out,outname,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
}

