library('bagfoot');
library('stringr')
library('argparse')

parser <- ArgumentParser(description='Export GRanges Object into Hotspot csv file like DNASE2Hotspots software ')
parser$add_argument("--bam",
                    type="character",
                    default=NULL,
                    help="bam file")

parser$add_argument("--outdir",
                    type="character",
                    default=NULL,
                    help="save result path")


parser$add_argument("--genome",
		    type="character",
                    default="mm10",
                    help="reference genome name")

parser$add_argument("--mappability",
		    type="character",
		    default=NULL,
		    help="MAPPABILITY_FILES_DIRECTORY")

args <- parser$parse_args()

############################### Configure
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}
ds239487_ptm <- proc.time();  # initial time
dlog <- function(output, ...) {
        etime = proc.time()-ds239487_ptm;
        cat(sprintf("[%g sec] %s", etime[["elapsed"]], output, ...));
        ds239487_ptm <- proc.time();
}


readBAMIndex<-function(bamFile, refgenome="mm9") {

	if (system("samtools", intern=F)!=1) {
		stop('Cannot run "samtools".  Install "samtools" first!');
	};

	if (refgenome=="mm9" | refgenome=="mm10" ) {
		chroms = paste("chr", c(1:19,"X","Y"), sep="");
	}

	if (refgenome =="hg19" | refgenome=="hg38") {
		chroms = paste("chr", c(1:22,"X","Y"), sep="");
	}

	baiFile = sprintf('%s.bai',bamFile);
	if (!isBigEnough(baiFile)) {
		dlog(sprintf("Generating the index file of %s\n", bamFile));
		system(sprintf('samtools index %s', bamFile));
	}
	idxout= system(sprintf('samtools idxstats %s', bamFile), intern=T);
	con <- textConnection(idxout);
	dat <- read.csv(con,sep='\t',header=F);
	names(dat) <- c("chr","maxloc","num","etc");
	dat[dat$chr %in% chroms,];
}

isBigEnough<-function(filename) {
	allowablesize = 300;
	if (file.exists(filename)) {
		ifelse(file.info(filename)$size > allowablesize , T,F);
	} else {
		F;
	}

}

MycountReadsBAM=function (bamFile, refgenome = "mm9")
{
    dat = readBAMIndex(bamFile, refgenome = refgenome)
    ncount = sum(dat$num)
    ncount
}



##############################
if(is.null(args$mappability)){
	MAPPABILITY_FILES_DIRECTORY=""
}else{
	MAPPABILITY_FILES_DIRECTORY=args$mappability 
}

##########  message dict ,will return 
message_dict=list()
##########  input paramaters
bamfile=args$bam
genome=args$genome

message_dict[["bamfile"]]=bamfile
message_dict[["genome"]]=genome
##########  count read bam
message("INFO : countReadsBAM ...")
cc = MycountReadsBAM(bamfile,genome);  # counts the number of cuts in the sequence file.
message_dict[["cc"]]=cc
########## makeCutCountBAM
message("INFO : makeCutCountBAM ...")
cutcountfile = makeCutCountBAM(bamfile,genome);   # generate a BedGraph file with DNase Cleavages counts
message_dict[["cutcountfile"]]=cutcountfile
########## MakeBiasCorrectionTableBAM
message("INFO : MakeBiasCorrectionTableBAM ...")
tabMappability= MakeBiasCorrectionTableBAM(bamfile= bamfile,
   outfile=paste0(cc,"_",genome,"_withMap.txt"),
   refgenome=genome,
   np=6,
   mapdir =MAPPABILITY_FILES_DIRECTORY);

message_dict[["tabMappability"]]=paste0(cc,"_",genome,"_withMap.txt")

