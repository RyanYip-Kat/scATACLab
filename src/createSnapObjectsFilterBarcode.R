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
parser$add_argument("--snapmeta",
                    type="character",
                    default=NULL,
                    help="meta csv file,include (sample_name,*.snp,barcode csv file)")



parser$add_argument("--outdir",
                    type="character",
                    default="Snap_result")

parser$add_argument("--addPmat",
		    action="store_true",
		    default=FALSE,
		    help="whether add Pmat")

args <- parser$parse_args()

#############################
outDir = args$outdir
metricDir=file.path(outDir,"metric")
mkdir(outDir)
mkdir(metricDir)


#############################
message("INFO : loading dataset ...")
DF=read.csv(args$snapmeta,sep=",",header=FALSE,stringsAsFactors=FALSE)
sample_list=as.character(DF$V1)
snap_list=as.character(DF$V2)
barcode_list=as.character(DF$V3)

message("INFO : create snp ...")
x.sp.list = lapply(seq_along(snap_list), function(i){ 
    cat(sprintf("INFO : create Snap : [ %d of %d ] --- [ %s ]\n",i,length(snap_list),sample_list[i]))
    x.sp = createSnap(file=snap_list[i], sample=sample_list[i],num.cores=8);
    ### select barcode
    barcodes = read.csv(barcode_list[i],head=TRUE,stringsAsFactors=FALSE)
    barcodes = barcodes[2:nrow(barcodes),];
    promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
    UMI = log(barcodes$passed_filters+1, 10);
    data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
    barcodes$promoter_ratio = promoter_ratio
    message("INFO : Plot metric ...")
    p1 = ggplot(
    		data, 
    		aes(x= UMI, y= promoter_ratio)) + 
                geom_point(size=0.1, col="grey") +
                theme_classic() +
                ggtitle(sample_list[i]) +
                ylim(0, 1) + xlim(0, 6) +
                labs(x = "log10(UMI)", y="promoter ratio") 
    ggsave(file.path(metricDir,paste0(sample_list[i],".pdf")),plot=p1,width=12,height=12)
        
    barcodes.sel = barcodes[which(UMI >= 3.2 & UMI <= 5 & promoter_ratio >= 0.25 & promoter_ratio <= 0.5),];
    barcodes.sel$barcode=paste(str_to_upper(sample_list[i]),as.character(barcodes.sel$barcode),sep="_")
    rownames(barcodes.sel) = barcodes.sel$barcode;
    x.sp = x.sp[which(x.sp@barcode%in%barcodes.sel$barcode),];
    x.sp@metaData = barcodes.sel[x.sp@barcode,]
    x.sp    
  })

names(x.sp.list) = sample_list

message("INFO : add bmat ...")
x.sp.list = lapply(seq_along(x.sp.list), function(i){
    cat(sprintf("INFO : addBmat : [ %d of %d ] --- [ %s ]\n",i,length(x.sp.list),sample_list[i]))
    x.sp = addBmatToSnap(x.sp.list[[i]], bin.size=5000);
    x.sp
  })

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

### add Pmat
if(args$addPmat){
	message("INFO : addPmat ...")
	x.sp = addPmatToSnap(x.sp);
}

x.sp@metaData$sample=x.sp@sample
message("INFO : Save ...")
saveRDS(x.sp,file.path(outDir,"snapObject.rds"))
message("INFO : Done!")


