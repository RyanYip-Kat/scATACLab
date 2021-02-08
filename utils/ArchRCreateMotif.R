library(ArchR)
library(Signac)
library(chromVAR)
library(motifmatchr)
library(chromVARmotifs)
library(SummarizedExperiment)
library(argparse)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
parser <- ArgumentParser(description='ChromVAR handle motif from ArchR')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="ArchR Project Path")


parser$add_argument("--outdir",
                    type="character",
                    default="output")
args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

outDir=args$outdir
makedir(outDir)
# getJasparMotifs() # "Mus musculus" or "Homo sapiens"
getMotifsDB=function (species = "Mus musculus", collection = "CORE", ...)
{
    opts <- list()
    opts["species"] <- species
    opts["collection"] <- collection
    opts <- c(opts, list(...))
    out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
    if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
        names(out) <- paste(names(out), TFBSTools::name(out),
            sep = "_")
    return(out)
}

message("INFO : Loading dataset ...")
projHeme=loadArchRProject(args$project)
features=getPeakSet(projHeme)
#MotifMatrix=getMatrixFromProject(projHeme,useMatrix="MotifMatrix")
#dev=new("chromVARDeviations", MotifMatrix)
#chromvar.z <- SummarizedExperiment::assays(dev)[[2]]
data("mouse_pwms_v2")
motifs = mouse_pwms_v2



MyAddMotifs <- function(
  object,
  genome,
  pfm,
  verbose = TRUE,
  ...
) {
  if (!requireNamespace("motifmatchr", quietly = TRUE)) {
    stop("Please install motifmatchr.\n",
         "https://www.bioconductor.org/packages/motifmatchr/")
  }
  if (verbose) {
    message("Building motif matrix")
  }
  motif.matrix <- CreateMotifMatrix(
    features = object,
    pwm = pfm,
    genome = genome,
    use.counts = FALSE
  )
  if (verbose) {
    message("Finding motif positions")
  }
  motif.positions <- motifmatchr::matchMotifs(
    pwms = pfm,
    subject = object,
    out = 'positions',
    genome = genome
  )
  if (verbose) {
    message("Creating Motif object")
  }
  motif <- CreateMotifObject(
    data = motif.matrix,
    positions = motif.positions,
    pwm = pfm
  )
  return(motif)
}


motif.matrix=MyAddMotifs(features,genome=BSgenome.Mmusculus.UCSC.mm10,pfm=mouse_pwms_v2)
#chromvar.z <- SummarizedExperiment::assays(dev)[[2]]
#rownames(x = chromvar.z) <- colnames(x = motif.matrix)
#obj <- CreateAssayObject(data = chromvar.z)
#MotifPlot(motif.matrix,motifs =c("Fosb","Mga"))
message("IFNO : Save...")
saveRDS(motif.matrix,file.path(outDir,"ArchRMotifObject.rds"))
message("INFO : Done!")
