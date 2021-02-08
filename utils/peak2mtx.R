library(argparse)
library(stringr)
library(Seurat)
library(DropletUtils)
library(Signac)
library(ArchR)

parser <- ArgumentParser(description='A Program to convert scATAC Peak matrix into  10X counts matrix...')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")


parser$add_argument("--object",
                    type="character",
                    default=NULL,
		    help="Signac,ArchRProject or SummarizedExperiment Object,rds file")


parser$add_argument("--assay",
		    type="character",
		    default="ATAC",
		    help="assay to use when object is Signac Object")
args <- parser$parse_args()

###############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

outDir=args$outdir
makedir(outDir)


message("INFO : Loading dataset ...")
proj=readRDS(args$object)

message("INFO : get Counts ...")
if(is(proj,"Seurat")){
	counts=GetAssayData(proj,slot="counts",assay=args$assay)
	peaks=granges(proj)
	peakDF=data.frame("chr"=seqnames(peaks),"start"=start(peaks),"end"=end(peaks))

}else if(is(proj,"SummarizedExperiment")){ # from ArchR
	counts=assay(proj)
	peaks=rowData(proj)$X
	peakDF=data.frame("chr"=seqnames(peaks),"start"=start(peaks),"end"=end(peaks))

}else if(is(proj,"ArchRProject")){
	peaks=getPeakSet(proj)
	se=getMatrixFromProject(proj,useMatrix="PeakMatrix")
	rowData(se)=peaks
	all.out=cbind(data.frame(chr=seqnames(peaks)),as.data.frame(ranges(peaks)))
	rownames(se)=peak
	counts=assay(se)

	peakDF=data.frame("chr"=seqnames(peaks),"start"=start(peaks),"end"=end(peaks))

}else{
	stop("INFO : Invalid Object!!!")
}


##################################
filename=file.path(outDir,"scATAC-count")
cat(sprintf("INFO : Write counts matrix into : [ %s ]\n",filename))
write10xCounts(x =counts, path=filename)
write.table(peakDF,file.path(filename,"peaks.bed"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
message("INFO : Done!")


