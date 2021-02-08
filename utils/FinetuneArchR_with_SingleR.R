library(argparse)
library(stringr)
library(SingleR)
library(ArchR)
library(Matrix)
library(ggplot2)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
		    default=NULL,
                    help="the path of project saved")

parser$add_argument("--singleR",
		    nargs="+",
                    type="character",
                    default="/home/ye/Data/dataset/SingleR-DataSet/SingleR-HumanPrimaryCellAtlasData.rds",
		    help="the singleR object as reference scRNA dataset")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--label",
                    type="character",
                    default="label.main",
                    help="the reference dataset'label (label.main,label.fine)")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)
print("# Get GeneScoreMatrix fro projHeme")

features=getFeatures(projHeme,useMatrix="GeneScoreMatrix")
meta=as.data.frame(getCellColData(projHeme))
experiment=getMatrixFromProject(projHeme,useMatrix="GeneScoreMatrix",verbose = TRUE)
rownames(experiment)=features

print("# Loading singleR dataset")
if(length(args$singleR)>1){
        ref_list<-lapply(args$singleR,function(ref){return(readRDS(ref))})
        label_list<-lapply(ref_list,function(ref){return(ref[[args$label]])})
}else{
        ref_list=readRDS(args$singleR)
        label_list=ref_list[[args$label]]
}

####################
safeBPParam <- function(nworkers) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}
(cores <- parallel::detectCores())
(bpParam <- safeBPParam(args$num_threads))


print("# Finetune with singleR")
mat=assay(experiment)
predictor <- SingleR(test = mat,
    ref = ref_list,
    labels =label_list,
    BPPARAM=bpParam)

print(table(predictor$labels))
pred.table=data.frame("Barcode"=colData(experiment)$Barcode,"singleR"=predictor$labels)
write.table(pred.table,file.path(args$outdir,"singleR-predictor.csv"),sep=",",row.names=F)
saveRDS(predictor,file.path(args$outdir,"singleR-predictor.rds"))
