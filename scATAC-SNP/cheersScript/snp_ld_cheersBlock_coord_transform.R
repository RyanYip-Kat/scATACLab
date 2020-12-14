library(argparse)
library(stringr)
library(liftOver)
parser <- ArgumentParser(description='coord transfrom,hg19 to hg38 or hg38 to hg19 or ...')
parser$add_argument("--chain",
                    type="character",
                    default="chainIO/hg19ToHg38.over.chain",
                    help="the over chain file")

parser$add_argument("--outdir",
                    type="character",
                    default="./results")

parser$add_argument("--snp",
                    type="character",
                    default="snp-list/Vogt-Koyanagi-Harada.txt",
                    help="SNP list file")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

message("INFO : Loading chain file")
chain <- import.chain(args$chain)
message("INFO : Loading SNP list file")
snp=read.table(args$snp,stringsAsFactors=FALSE,sep="\t",header=TRUE)  # Chrom,SNP,BP
snp$start=snp$BP-1

message("INFO : make GRanges")
gr=makeGRangesFromDataFrame(snp,ignore.strand=T,seqnames.field="Chrom",end.field="BP",start.field="start")
mcols(gr)$rsid=snp$SNP
message("INFO : transform coord")
gr_transform <- liftOver(gr,chain)

message("INFO : save new snp after coord transform")
new_snp=as.data.frame(gr_transform)
new_snp=new_snp[,c("seqnames","rsid","end")]
colnames(new_snp)=c("Chrom","SNP","BP")
new_snp=na.omit(new_snp)
write.table(new_snp,file.path(args$outdir,basename(args$snp)),sep="\t",row.names=FALSE,quote=FALSE)



