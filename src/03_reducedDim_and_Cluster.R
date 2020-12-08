library(Seurat)
library(Signac)
library(SeuratWrappers)
library(harmony)
library(argparse)
library(stringr)
library(future)
library(ggplot2)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")


parser$add_argument("--outdir",
                    type="character",
                    default="./Results")

parser$add_argument("--batch_key",
                    type="character",
                    default=NULL,
                    help="key used for batch correct")

parser$add_argument("--batch_correct",
                    action='store_true', default=FALSE)


parser$add_argument("--method",
                    type="character",
                    default="harmony",
                    choices=c("fastmnn","harmony"),
                    help="batch correct method")


args <- parser$parse_args()

############################### Configure
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

plan("multiprocess", workers = 16)

outDir=args$outdir
makedir(outDir)
############################## 

message("INFO : Loading dataset")
object=readRDS(args$seurat)
##############################
message("INFO : Run LIS")
object <- RunTFIDF(object)
object <- FindTopFeatures(object, min.cutoff = 'q0')
object <- RunSVD(object)


############################### batch correct
use_rep="lsi"
dims.use=2:30
if(args$batch_correct){
        vars=ifelse(!is.null(args$batch_key),args$batch_key,"Sample")
        stopifnot(vars%in%colnames(object@meta.data))
        if(args$method=="harmony"){
                cat(sprintf("INFO : Run Harmony with [ %s ] \n",vars))
                object=RunHarmony(object, group.by.vars =vars,
                                  reduction ="lsi",
				  dims.use=dims.use,
				  max.iter.harmony=20,
				  assay.use=DefaultAssay(object))
		use_rep="harmony"
                dims.use=1:30

        }else if(args$method=="fastmnn"){
		message("INFO : Run fastmnn method")
                object.list = SplitObject(object, split.by = vars)
                object=RunFastMNN(object.list,features=10000,assay=DefaultAssay(object))
                use_rep="mnn"
		dims.use=1:30
        }else{
                stop("Invalid batch correct method !")
        }
}

#################################  Run Reduce Embeddings
message("INFO : Run UMAP ...")
object <- RunUMAP(object = object, reduction =use_rep, dims = dims.use)

message("INFO : Clustering ...")
object <- FindNeighbors(object = object, reduction =use_rep, dims = dims.use)
object <- FindClusters(object = object, verbose = FALSE, algorithm = 3)

#################################
message("INFO : Save Object ...")
saveRDS(object,file.path(outDir,"seurat.rds"))


################################# Dim Plots
plotDir=file.path(outDir,"Plots")
makedir(plotDir)
p=DimPlot(object = object, label = TRUE) + NoLegend()
ggsave(file.path(plotDir,"Clusters.pdf"),plot=p,width=12,height=10)
message("INFO : Done!")
