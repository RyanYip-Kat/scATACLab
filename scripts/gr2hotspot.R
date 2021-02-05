library(ArchR)
library(stringr)
library(argparse)

parser <- ArgumentParser(description='Export GRanges Object into Hotspot csv file like DNASE2Hotspots software ')
parser$add_argument("gr",
                    type="character",
                    default=NULL,
                    help="Granges Object rds file")

parser$add_argument("outdir",
                    type="character",
                    default=NULL,
                    help="save result path")

args <- parser$parse_args()

############################### Configure
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

write2hotspot<-function(file,outdir,name="hotspot"){
	gr=readRDS(file)
	chr=as.character(seqnames(gr))
	st=as.integer(start(gr))
	ed=as.integer(end(gr))
	nr=length(chr)
	hotspot = data.frame(ID=1:nr, chrom.=chr, start= st + 1, end = ed,
                        MaxD = rep(0, nr),
                        AveD = rep(0, nr),
                        Zscore = rep(0, nr),
                        pvalue = rep(0, nr));

	makedir(outdir)
	if(is.null(name)){
		name=str_split(basename(file),".rds")[[1]][1]
	}
	csvfile= file.path(outdir,paste0(name,".csv"));
        write.table(hotspot, file=csvfile,sep=",", row.names=F);
        cat(sprintf('INFO : writing  --- %s\n',csvfile));
}
##############################
outDir=args$outdir
grFile=args$gr
message("INFO : Get hotspot file ...")
write2hotspot(file=grFile,outdir=outDir,name=NULL)
message("INFO : Done!")


