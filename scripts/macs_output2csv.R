#!/usr/bin/Rscript
# written by Songjoon BAek
# 09.08.2016 10:44:16

args <- commandArgs(trailingOnly = TRUE);

#args = '/home/baeks/project/gfoot/metab_new/macs/output_fasted/fasted_merged_peaks.narrowPeak';

options("browserNLdisabled"=T);

#source("~/Data/R/tools/HotspotUtil.R")

if (length(args)<1) {
	stop("Usage:  macs_output2csv [MACS2 narrowPeak file]\n");
}


hasSubstring<- function(string, pattern)  {
	result= length(grep(tolower(pattern), tolower(string)))>0;
	result;
}

#
infiles= Sys.glob(args[1]);

#infile = infiles[1]
for (infile in infiles) {

	if (hasSubstring(infile, ".narrowPeak")) {
		macs = read.csv(infile,sep="\t",header=F,stringsAsFactors=F);     # hotspot csv file name.
		names(macs)[1:3] <- c('chr','st','ed');
		nr = nrow(macs);
		
		hotspot = data.frame(ID=1:nr, chrom.=macs$chr, start= macs$st + 1, end = macs$ed, 
			MaxD = rep(0, nr),
			AveD = rep(0, nr),
			Zscore = rep(0, nr),
			pvalue = rep(0, nr));
				   
		csvfile= file.path(getwd(),gsub('.narrowPeak','.csv',basename(infile)));
		write.table(hotspot, file=csvfile,sep=",", row.names=F);
		cat(sprintf('writing %s\n',	csvfile));
	} else {
		stop('Not supported file.');
	}
}
