library(gkmSVM) 
library(argparse)
library(stringr)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--path",
                    type="character",
                    default=NULL,
                    help="The bed file for generating pos,neg fa file")

parser$add_argument("--outdir",
                    type="character",
                    default="./result")

args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



bedFiles=list.files(args$path,"*.bed")
message("### Generate Seqs Response to Peak")
for(bed in bedFiles){
	name=str_split(bed,"-")[[1]][1]
	bedfile=file.path(args$path,bed)
	outfile=file.path(args$outdir,name)

	if(!dir.exists(outfile)){
		dir.create(outfile,recursive=TRUE)
	}
	cat(sprintf("### Generate sequence for %s --- from : %s\n",name,bedfile))
	genNullSeqs(bedfile,nMaxTrials=10,
	    xfold=3,
	    genomeVersion='hg38',
	    outputPosFastaFN=file.path(outfile,'ctcfpos.fa'),
	    outputBedFN=file.path(outfile,'ctcf.bed'), 
	    outputNegFastaFN=file.path(outfile,'ctcfneg.fa'))
}

