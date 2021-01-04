library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(yaml)
library(Rcpp)
library(argparse)
library(SnapATAC)
library(stringr)

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


# gtfFile <- "data/genes.gtf"
#genes <- getGeneGTF(gtfFile) %>% resize(1,"start") %>% resize(tssWindow * 2 + 1, "center")
##############################
parser <- ArgumentParser(description='Create snap Object from mutiple sample ...')
parser$add_argument("--snap",
                    type="character",
                    default=NULL,
                    help="snap object file from create Object")



parser$add_argument("--outdir",
                    type="character",
                    default="Snap_result")


parser$add_argument("--bed",
                    type="character",
                    default="mm10_Gencode_VM18.bed.gz",
                    help="gencode bed file file")

args <- parser$parse_args()

#############################
outDir = args$outdir
bedFile=args$bed
mkdir(outDir)


############################
message("INFO : loading dataset ...")
x.sp=readRDS(args$snap)

message("INFO : createGmatFromMat ...")
genes.df = read.table(bedFile);
genes.gr = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4])
x.sp = createGmatFromMat(
    obj=x.sp,
    input.mat="bmat",
    genes=genes.gr,
    do.par=TRUE,
    num.cores=10
);

message("INFO : Save ...")
saveRDS(x.sp,file.path(outDir,"snapAnnotation.rds"))
message("INFO : Done!")
