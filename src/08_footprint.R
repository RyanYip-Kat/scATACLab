library(argparse)
library(motifmatchr)
library(stringr)
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
outDir=args$outdir
genome=args$genome
makedir(outDir)


############################## Loading dataset
message("INFO : Loading dataset ...")
seurat=readRDS(args$seurat)


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

##############################
message("INFO : extract position frequency matrices for the motifs ...")
pwm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = species, all_versions = FALSE)
)



message("INFO : match motif positions ...")
motif.positions <- matchMotifs(
  pwms = pwm,
  subject = granges(seurat),
  out = 'positions',
  genome = genome
)


message("INFO : create a Motif object and add it to the assay ...")
motif <- CreateMotifObject(
  positions = motif.positions,
  pwm = pwm
)

seurat <- SetAssayData(
  object = seurat,
  slot = 'motifs',
  new.data = motif
)

message("INFO : Save Object ...")
saveRDS(seurat,file.path(outDir,"seurat-footprint.rds"))
message("INFO : Done!")
