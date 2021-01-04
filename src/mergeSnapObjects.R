library(SnapATAC)
library(GenomicRanges)
library(viridisLite)
library(ggplot2)
library(argparse)
library(stringr)
############################
mkdir=function(path){
	if(!dir.exists(path)){
		dir.create(path,recursive=TRUE)
	}
}


############################
parser <- ArgumentParser(description='Create snap Object from mutiple sample ...')
parser$add_argument("--snaps",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="snapObject.rds files")

parser$add_argument("--outdir",
                    type="character",
                    default="Snap_result")

args <- parser$parse_args()

#############################
outDir = args$outdir
snap_files=args$snaps
mkdir(outDir)


#############################
message("INFO : loading dataset ...")
x.sp.list=lapply(seq_along(snap_files),function(i){
			 file=snap_files[i]
			 cat(sprintf("INFO : loading [ %d of %d ] --- [ %s ]\n",i,length(snap_files),file))
			 x.sp=readRDS(file)
			 return(x.sp)})



message("INFO : combine snap objects ...")
bin.shared = Reduce(intersect, lapply(x.sp.list, function(x.sp) x.sp@feature$name));
x.sp.list <- lapply(x.sp.list, function(x.sp){
    idy = match(bin.shared, x.sp@feature$name);
    x.sp[,idy, mat="bmat"];
  })
x.sp = Reduce(snapRbind, x.sp.list);

message("INFO : Binarize matrix ...")
x.sp = makeBinary(x.sp, mat="bmat");

message("INFO : Filter bins ...")
### filter bin
black_list = read.table("mm10.blacklist.bed.gz");
black_list.gr = GRanges(
    black_list[,1],
    IRanges(black_list[,2], black_list[,3])
  );
idy = queryHits(findOverlaps(x.sp@feature, black_list.gr));
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};

### remove unwanted chroms
chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))];
idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature);
if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};

### the bin coverage roughly obeys a log normal distribution. We remove the top 5% bins that overlap with invariant features such as promoters of the house keeping genes
bin.cov = log10(Matrix::colSums(x.sp@bmat)+1);
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95);
idy = which(bin.cov <= bin.cutoff & bin.cov > 0);
x.sp = x.sp[, idy, mat="bmat"];

message("INFO : Save ...")
saveRDS(x.sp,file.path(outDir,"snapObject.rds"))
message("INFO : Done!")


