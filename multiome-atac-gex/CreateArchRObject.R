library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--inputFiles",
                    type="character",
                    default=NULL)


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--genome",
                    type="character",
                    default="mm10",
		    choices=c("hg38","hg19","mm9","mm10"),
		    help="which genome to be used")

parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")


parser$add_argument("--rm_doublet",
		    action="store_true",
                    default=FALSE,
                    help="whether to remove doublet")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


inputFiles<-read.table(args$inputFiles,sep=",",stringsAsFactors=FALSE)
files=as.character(inputFiles$V1)
h5Files=as.character(inputFiles$V2)
SampleNames=as.character(inputFiles$V3)
Tss=as.numeric(inputFiles$V4)
Frags=as.numeric(inputFiles$V5)

print(paste0("Setting default genome  with :",args$genome))
addArchRGenome(args$genome)

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 


###################
max_frags=50000
max_tss=15

#################
print("# Create ArrowFiles")
ArrowFiles=c()
for(i in 1:length(files)){
	arrow=createArrowFiles(
			       inputFiles = files[i],
                               sampleNames =SampleNames[i],
                               filterTSS = Tss[i], #Dont set this too high because you can always increase later
                               filterFrags =Frags[i],
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
			       excludeChr = c("chrM", "chrY"),
			       minFrags = 500,
			       maxFrags = max_frags,
			      )
	ArrowFiles=c(ArrowFiles,arrow)
	print(arrow)
}
print(ArrowFiles)

message("INFO : import 10x FeatureMatrix ...")
seRNA <- import10xFeatureMatrix(
    input = h5Files,
    names = SampleNames
)
if(length(seRNA)>1){
	seRNA=do.call(cbind,seRNA)
}

saveRDS(seRNA,file.path(args$outdir,"seRNA.rds"))
message("INFO : Create Object ...")
projHeme <- ArchRProject(
  ArrowFiles = ArrowFiles,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

cellPass=intersect(colnames(seRNA),getCellNames(projHeme))
idx=getCellNames(projHeme)%in%cellPass
projHeme=projHeme[idx]
idx=colnames(seRNA)%in%cellPass
seRNA=seRNA[,idx]

message("INFO : add scRNA ...")
projHeme <- addGeneExpressionMatrix(input =projHeme, seRNA = seRNA, force = TRUE,excludeChr = c("chrM", "chrY"))

#if(is.null(getImputeWeights(projHeme))){
#        projHeme=addImputeWeights(projHeme)
#}

print("### Save ArchR Object")
saveArchRProject(ArchRProj = projHeme, outputDirectory =file.path(args$outdir,"Save-ProjHeme-Create"), load =TRUE)
#saveRDS(projHeme,file.path(args$outdir,"Save-ProjHeme-Create/Save-ArchR-Project.rds")
#subsetArchRProject(projHeme,cells=cellPass,outputDirectory=file.path(args$outdir,"Save-ProjHeme-Create2"))
if(args$rm_doublet){
    message("INFO : get doublet score ...")
    cmd=paste0("bash /home/ye/Work/R/scATAC/ArchR/multiome-atac-gex/getDoubletScore.sh"," ", file.path(args$outdir,"Save-ProjHeme-Create")," ",args$outdir)
    system(cmd)

    #DF=read.csv(file.path(args$outdir,"DoubletScores.csv"))
    #for(name in colnames(DF)){
    #	cat(sprintf("INFO : %s \n",name))
    #	projHeme=addCellColData(projHeme,data=DF[,name],name=name,cells=rownames(DF))
    #}
    #projHeme=filterDoublets(projHeme)
    #saveRDS(projHeme,file.path(args$outdir,"Save-ProjHeme-Create/Save-ArchR-Project.rds"))
}


