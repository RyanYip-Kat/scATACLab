library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
library(ggplot2)
library(edgeR)
library(matrixStats)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--se",
                    type="character",
                    default=NULL,
		    help="The object from getMarkerFeatures")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--fdr",
		    type="character",
                    default="0.1")

parser$add_argument("--fc",
                    type="character",
                    default="0.25")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

se=readRDS(args$se)
cutOff=paste0("FDR < ",args$fdr," & ","Log2FC > ",args$fc)
cat(sprintf("### The cutoff condition is :%s\n",cutOff))
markerList <- getMarkers(se, cutOff = cutOff,returnGR=TRUE)
for(name in names(markerList)){
	path=file.path(args$outdir, paste0(str_replace(name," ","_"), "-Peak.bed"))
	fragmentsj <- markerList[[name]]
	cat(sprintf("### Write : %s Different peaks into : %s \n",name,path))
	#message(sprintf("### Write : %s Different peaks into : %s bed file\n",name,path))
	if(length(fragmentsj) > 0){
	    #out <- data.frame(
	    #  chr = c(seqnames(fragmentsj), seqnames(fragmentsj)),
	    #  start = c(as.integer(start(fragmentsj) - 1), as.integer(end(fragmentsj) - 1)),
	    #  end = c(as.integer(start(fragmentsj)), as.integer(end(fragmentsj)))
	    #  ) 
	  out=cbind(data.frame(chr=seqnames(fragmentsj)),
		    as.data.frame(ranges(fragmentsj))
		    )
	  out=out[,c("chr","start","end")]
	  readr::write_tsv(
		    x = out,
                    append = FALSE,
                    path = path,
                    col_names = FALSE)
	}

}
