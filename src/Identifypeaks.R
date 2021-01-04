library(Matrix)
library(GenomicRanges)
library(SummarizedExperiment)
library(argparse)
library(SnapATAC)
library(stringr)

###############################
mkdir=function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

findMacs2=function(){
        macs2_path=system("which macs2",intern=TRUE)
        return(macs2_path)
}

findSnaptools=function(){
        path=system("which snaptools",intern=TRUE)
        return(path)
}

##############################
parser <- ArgumentParser(description='Create snap Object from mutiple sample ...')
parser$add_argument("--snap",
                    type="character",
                    default=NULL,
                    help="snap object file from create Object")

parser$add_argument("--callCol",
		    type="character",
		    default="cluster",
		    help="which column for callpeak")

parser$add_argument("--gsize",
		    type="character",
		    default="mm",
		    choices=c("mm", "hs"),
		    help="gsize for call peak")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
                    help="column in snap metaData")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL)


parser$add_argument("--outdir",
                    type="character",
                    default="Snap_result")


args <- parser$parse_args()

#############################
outDir = args$outdir
Temp=file.path(outDir,"tmp")
mkdir(outDir)
mkdir(Temp)


#############################
message("INFO : loading dataset ...")
x.sp=readRDS(args$snap)
metadata=x.sp@metaData
if(!is.null(args$column) & !is.null(args$subset)){
        assertthat::assert_that(args$column%in%colnames(metadata))
        metadata[[args$column]]=as.character(metadata[[args$column]])
        column=metadata[[args$column]]
        idx=args$subset%in%unique(column)
        assertthat::assert_that(sum(idx)==length(args$subset))
        idx=which(column%in%args$subset)
        x.sp=x.sp[idx,]
}

#############################
message("INFO : call peak ...")
metadata=x.sp@metaData
# call peaks for all cluster with more than 100 cells
clusters.sel = names(table(metadata[[args$callCol]]))[which(table(metadata[[args$callCol]]) > 100)];
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i]);
    peaks = runMACS(
        obj=x.sp[which(metadata[[args$callCol]]==clusters.sel[i]),], 
        output.prefix=paste0("Snap.", gsub(" ", "_", clusters.sel)[i]),
        path.to.snaptools=findSnaptools(),
        path.to.macs=findMacs2(),
        gsize=args$gsize, # mm, hs, etc
        buffer.size=500, 
        num.cores=1,
        macs.options="--nomodel --shift 100 --ext 200 --qval 5e-2 -B --SPMR",
        tmp.folder=Temp
   );
	peaks
 	}, mc.cores=12);


##############################
message("INFO  : combined peaks ...")
peaks.names = system("ls | grep narrowPeak", intern=TRUE);
peak.gr.ls = lapply(peaks.names, function(x){
    peak.df = read.table(x)
    GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
  })
peak.gr = reduce(Reduce(c, peak.gr.ls));
peaks.df = as.data.frame(peak.gr)[,1:3];

##############################
message("INFO : Save ...")
write.table(peaks.df,file = file.path(outDir,"peaks.combined.bed"),append=FALSE,
		quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".", 
		row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
		fileEncoding = "")

message("INFO : Done!")
