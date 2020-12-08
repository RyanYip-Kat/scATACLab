library(argparse)
library(stringr)
library(future)
library(Signac)
library(Seurat)
library(JASPAR2018)
library(TFBSTools)
library(patchwork)
library(ggplot2)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")

parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")


parser$add_argument("--outdir",
                    type="character",
                    default="./Results")


args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

findMacs2=function(){
	macs2_path=system("which macs2",intern=TRUE)
	return(macs2_path)
}

############################## Configure
plan("multiprocess", workers = 16)
outDir=args$outdir
genome=args$genome
makedir(outDir)


############################## Loading dataset
message("INFO : Loading dataset")
object=readRDS(args$seurat)


##############################
if(genome=="hg38"){
	library(BSgenome.Hsapiens.UCSC.hg38)
	BS=BSgenome.Hsapiens.UCSC.hg38
	species="Homo sapiens"
}else{
	library(BSgenome.Mmusculus.UCSC.mm10)
	BS=BSgenome.Mmusculus.UCSC.mm10
	species="Mus musculus"
}


##############################  add Motif
# Get a list of motif position frequency matrices from the JASPAR database
message("INFO : getMatrixSet ...")
pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species =species, all_versions = FALSE)
)

# add motif information
message("INFO : Add Motifs ...")
object <- AddMotifs(
  object = object,
  genome = BS,
  pfm = pfm
)


############################### Computing motif activities
message("INFO : Computing motif activities")
object <- RunChromVAR(
  object = object,
  genome = BS
)
message("INFO : Save Object")
saveRDS(object,file.path(outDir,"seurat-motif.rds"))
message("INFO : Done!")

p=MotifPlot(object,motifs="MA0497.1",assay = 'peaks')
ggsave(file.path(outDir,"MA0497.1.pdf"),plot=p,width=12,height=10)

