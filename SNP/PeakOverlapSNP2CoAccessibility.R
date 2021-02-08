library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
library(ggplot2)
library(edgeR)
library(matrixStats)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the project path  of ArchR")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--binarize",
                    action="store_true")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
		    help="which column to subset")

parser$add_argument("--subset",
                    type="character",
		    nargs="+",
                    default=NULL)

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
corCutOff=0.1
projHeme=loadArchRProject(args$project)
metadata=as.data.frame(getCellColData(projHeme))
stopifnot(args$groupby%in%colnames(metadata))
stopifnot(args$useMatrix%in%getAvailableMatrices(projHeme))

cells=getCellNames(projHeme)
if(!is.null(args$column) & !is.null(args$subset)){
	stopifnot(args$column%in%colnames(metadata))
	target=metadata[[args$column]]
	
	stopifnot(args$subset%in%unique(target))
	idxPass= which(target%in%args$subset)
	cellPass=cells[idxPass]
	projHeme=subsetCells(projHeme,cellNames=cellPass)
}
###############
message("### Get CoAccessibility")
cA=getCoAccessibility(projHeme,corCutOff=corCutOff,returnLoops=TRUE)  #   co-accessible above 0.3
#projHeme5 <- addPeakAnnotations(ArchRProj = projHeme5, regions =c("SNP"="snp.bed"), name = "SNP",force=T)  # it worked

snp_list="/home/ye/Work/BioAligment/SNP/Shi/masterfile9.csv"
df=read.csv(snp_list,stringsAsFactors=F,header=TRUE)
gr=makeGRangesFromDataFrame(df,seqnames.field=c("seqnames", "seqname","chromosome", "chrom","chr", "chromosome_name","seqid","Chr"),
						,start.field="pos",end.field="pos")

mcols(gr)=df
diseases=unique(mcols(gr)$Disease)
cols=c("chr","pos","pos")
beds=c()
for(disease in diseases){
	df_subset=subset(df,Disease==disease)
	df_subset=df_subset[,cols]
	write.table(df_subset,file.path(args$outdir,paste0(disease,".bed")),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	x=file.path(args$outdir,paste0(disease,".bed"))
	names(x)=disease
	beds=c(beds,x)
}

projHeme <- addPeakAnnotations(ArchRProj = projHeme, regions = beds, name = "SNP")
projHeme=addDeviationsMatrix(projHeme,peakAnnotation="SNP")

regionPositions=GRangesList(lapply(diseases,function(x)gr_disease=subset(gr,Disease==x)))
names(regionPositions) <- diseases
peakSet <- getPeakSet(projHeme)
allPositions <- unlist(regionPositions)

regionMat <- Matrix::sparseMatrix(
    i = queryHits(overlapRegions),
    j = match(names(allPositions),names(regionPositions))[subjectHits(overlapRegions)],
    x = rep(TRUE, length(overlapRegions)),
    dims = c(length(peakSet), length(regionPositions))
  )
colnames(regionMat) <- names(regionPositions)
regionMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = regionMat,counts=regionMat), rowRanges = peakSet)
out <- SimpleList(
      regionMatches = regionMat,
      regionPositions = regionPositions,
      date = Sys.Date()
    )

message("INFO : Save Results")
saveRDS(out,file.path(args$outdir,"regionMatches.rds"))
snpMatrix=getMatrixFromProject(proj,useMatrix="SNPMatrix")
saveRDS(snpMatrix, file.path(args$outdir,"GWAS-addSNPMatrix-chromVAR-Summarized-Experiment.rds"))
#saveArchRProject(ArchRProj = projHeme, outputDirectory = file.path(args$outdir,"Save-ProjHeme-SNPMatrix"), load = TRUE)
