library('bagfoot');

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

calcFrequencyTableBAMChromosome<-function(bamfile, chr, genome=NA, np=6,mapdir='',shifts=c(0,0)) {

  tempdata = sprintf('temp_data_%s.dat',digest::digest(list(bamfile, chr, GenomeInfoDb::providerVersion(genome),np, shifts, mapdir)));
  if (file.exists(tempdata)) {
#  	  cat(sprintf('calcFrequencyTableBAMChromosome: loading saved data for %s: refgenome=%s\n', chr, providerVersion(genome)));
    freqtable = readRDS(tempdata);
  } else {

#~ 	  if (mappability) {
#~ 		 mapdir = getMapDirectory(genome);
#~ 	  } else {
#~ 		  mapdir ="";
#~ 	  }

	  #unmappable = pickUnmappaleBasesByMappability(mapdir, chr);

	  cuts=readCutSitesPerChromFromBAM(bamfile, chr, shifts=shifts);  # ok

	 # cat(sprintf('calcFrequencyTableBAMChromosome: finished reading cuts for %s\n', chr));
	  nuccode<-readNucleotideCodeForChromosome(chr,nmer=np,genome=genome,mapdir=mapdir); # ok


	  cutfreq= rle(cuts);
	  tab=cbind(loc=cutfreq$values, count = cutfreq$lengths);
	  nucfreq = vector("numeric", length=length(nuccode$freqtable$count));
	  codetochange=nuccode$code[cutfreq$values-floor(np/2)];
	  for (x in 1:length(codetochange)) {
		nucfreq[codetochange[x]] = nucfreq[codetochange[x]] + cutfreq$lengths[x];
	  }

      freqtable = data.frame(count=nucfreq,refseqcount= nuccode$freqtable$count, row.names=nuccode$freqtable$seq);
      cat(sprintf('calcFrequencyTableBAMChromosome: saving the result for %s\n', chr));
	  saveRDS(freqtable, file = tempdata);
	 # browser()
  }
  freqtable;
}

calcFreqeuncyTableBAM <- function ( bamfile='',refgenome='', np=6,mapdir='',shifts=shifts) {

  tempdata = sprintf('temp_main_data_%s.dat',digest::digest(list(bamfile, refgenome,np, shifts, mapdir)));
  if (file.exists(tempdata)) {
    freqtable = readRDS(tempdata);
  } else {

    genome <- NA;
    chrominfo = loadChromosomeRange(refgenome);
    chroms = names(chrominfo);
    genome <- loadReferenceGenome(refgenome);

   # print(chroms);

   freqs=parallel::mclapply(chroms, function(x) calcFrequencyTableBAMChromosome(bamfile, x,genome=genome, np=np,mapdir=mapdir,shifts=shifts), mc.cores=10);
  #freqs=lapply(c('chr1','chr2'), function(x) calcFrequencyTableBAMChromosome(bamfile, x,genome=genome, np=np,mapdir=mapdir,shifts=shifts));

    countsum = freqs[[1]]$count;
    refseqcountsum = freqs[[1]]$refseqcount;
    if (length(freqs)>1) {
      for (ll in 2:length(freqs)) {
        countsum = countsum + freqs[[ll]]$count;
        refseqcountsum = refseqcountsum + freqs[[ll]]$refseqcount;
      }
    }

    freqtable = data.frame(count=countsum,refseqcount= refseqcountsum, row.names=row.names(freqs[[1]]));
    saveRDS(freqtable, file = tempdata);
  }
  freqtable;
}

loadReferenceGenome<- function(refgenome) {
	genome = NULL;
	if (refgenome=='mm9') {
      genome <- BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9;
    } else if (refgenome=='hg19') {
      genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19;
    } else if (refgenome=='mm10') {
	  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10;
	} else if (refgenome=='hg38') {
	  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38;
	} else if (refgenome=='rn5') {
	  genome <- BSgenome.Rnorvegicus.UCSC.rn5::BSgenome.Rnorvegicus.UCSC.rn5;
	} else if (refgenome=='rn6') {
	  genome <- BSgenome.Rnorvegicus.UCSC.rn6::BSgenome.Rnorvegicus.UCSC.rn6;
	} else {
	   cat('Supported Ref. genomes : mm9, mm10, hg19, hg38, rn5, rn6\n');
	   stop(sprintf('not supported ref. genome %s\n', refgenome));
	}
    genome;
}

loadChromosomeRange<-function(refgenome) {
# last edited: 2/2/2018

  allowableChrs = c(paste('chr',1:22,sep=''),"chrX","chrY");
  genome=loadReferenceGenome(refgenome);
  seqnames = genome@seqinfo@seqnames;
  filter= seqnames %in% allowableChrs;
  seqnames = seqnames[filter];
  seqlengths= genome@seqinfo@seqlengths[filter];
  ll <- lapply(seq_along(seqnames), function(x) c(as.integer(1), as.integer(seqlengths[x])));
  names(ll) = seqnames;
  ll;
}

##############################
MAPPABILITY_FILES_DIRECTORY_MM10="mm10/50mers"
bamfile1="/home/ye/Data/dataset/cellranger-aggr/CEpi-ATAC/mouse/MSCE-24h-1-ATAC/outs/possorted_bam.bam"
cc1 = MycountReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.

cutcountfile1 = makeCutCountBAM(bamfile1,"mm10");   # generate a BedGraph file with DNase Cleavages counts
tabMappability10 = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="cepi_24h1_mm10_withMap.txt",
   refgenome="mm10",
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM10);
