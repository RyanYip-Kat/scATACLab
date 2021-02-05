library(ArchR)
library(stringr)
library(argparse)

parser <- ArgumentParser(description='Export GRanges Object into Hotspot csv file like DNASE2Hotspots software ')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="ArchR Object rds file")

parser$add_argument("--outdir",
                    type="character",
                    default=NULL,
                    help="save result path")


parser$add_argument("--split",
		    action="store_true",
		    default=FALSE,
		    help="whether split Granges Object from ArchR and export them respect")

args <- parser$parse_args()

############################### Configure
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

write2hotspot<-function(gr,outdir,name=NULL){
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
		name="hot2pot"
	}
	csvfile= file.path(outdir,paste0(name,".csv"));
        write.table(hotspot, file=csvfile,sep=",", row.names=F);
        cat(sprintf('INFO : writing  --- %s\n',csvfile));
}
##############################
outDir=args$outdir
project=args$project
name=args$name
split=args$split

message("INFO : Get GRanges Object ...")
projHeme=loadArchRProject(project,showLogo=FALSE)
GR=getPeakSet(projHeme)
message("INFO : Get hotspot file ...")
if(split){
	GRNames=unique(names(GR))
	for(i in seq_along(GRNames)){
		grName=GRNames[i]
		cat(sprintf("INFO : [ %d ( %s ) of %d ]\n",i,grName,length(GRNames)))
		gr=GR[names(GR)%in%grName]
		write2hotspot(gr=gr,outdir=outDir,name=grName)
	}
}else{
        write2hotspot(gr=GR,outdir=outDir,name="ArchRPeakSet")
}
message("INFO : Done!")


