#Running gwas chromVAR on single-cell summarized experiment using co-accessibility
#05/02/19
#Cite Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(ArchR)
library(chromVAR)
library(Matrix)
#library(MatrixStats)
library(SummarizedExperiment)
library(magrittr)
library(BiocParallel)
library(argparse)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
register(SerialParam())
set.seed(1)

parser <- ArgumentParser(description='Running gwas chromVAR on single-cell summarized experiment using co-accessibility')
#parser$add_argument("--project",
#                    type="character",
#                    default=NULL,
#                    help="the project path  of ArchR")
parser$add_argument("--se",
		    type="character",
		    default=NULL,
		    help="PeakSet SummarizedExperiment Object")

parser$add_argument("--outdir",
                    type="character",
                    default="./results")


#parser$add_argument("--binarize",
#                    action="store_true")

parser$add_argument("--conn",
                    type="character",
                    default=NULL,
		    help="Peaks-Co-Accessibility from cicero")


parser$add_argument("--cutoff",
		    type="double",
		    default=0.35,
		    help="The cutoff value for Co-Accessibility filtering")

parser$add_argument("--SNP",
                    type="character",
                    default=NULL,
                    help="GWAS SNP GRange object")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

message("INFO : Get Peak Set ")
se=readRDS(args$se)
#projHeme=loadArchRProject(args$project)
#se=getMatrixFromProject(projHeme,useMatrix="PeakMatrix")
#peaks=getPeakSet(projHeme)
#rowRanges(se)=peaks
#rowData(se)=cbind(rowData(se),mcols(peaks))
#rgs=as.data.frame(ranges(peaks))
#sqs=as.data.frame(seqnames(peaks))
#peakset=cbind(sqs,rgs)
#peak_set=paste(peakset$value,peakset$start,peakset$end,sep="_")
#rownames(se)=peak_set
#print(head(rowRanges(se)))

message("INFO : Loading GWAS SNP and Co-Accessibility")
gr <- readRDS(args$SNP)
conn <- readRDS(args$conn)
conn <- conn[conn[,3] >= args$cutoff,]
peaknames <- paste(seqnames(se),start(se),end(se),sep="_")
conn[,4] <- match(paste0(conn[,1]), peaknames)
conn[,5] <- match(paste0(conn[,2]), peaknames)
connMat <- Matrix::sparseMatrix(i=conn[,4],j=conn[,5],x=rep(1,nrow(conn)),dims=c(nrow(se),nrow(se)))

#Add Bias
message("INFO : Add Bias")
genome <- BSgenome.Hsapiens.UCSC.hg38
se <- addGCBias(se, genome = genome)

#Overlap GWAS SNPs
message("INFO Overlap GWAS SNPs")
o <- lapply(split(gr, gr$Disease), function(x){
	#extend snp +- 500 bp
	overlapsAny(se, resize(x,1001,"center"), ignore.strand = TRUE)
}) %>% Reduce("cbind",.)
ow <- which(o > 0, arr.ind=TRUE)
matches <- Matrix::sparseMatrix(i=ow[,1],j=ow[,2],x=o[cbind(ow[,1],ow[,2])], dims = c(nrow(se),length(unique(gr$Disease))))
colnames(matches) <- names(split(gr, gr$Disease))

#Use connections mat!
message("INFO : Use connections match")
matches2 <- matches
idx <- which(Matrix::rowSums(matches2)>0) #which peaks have a snp
for(i in seq_along(idx)){
	if(i %% 100 == 0){message(sprintf("%s of %s",i,length(idx)))}
	#peaks co-accessible to peak with snp
	coi <- unique(c(which(connMat[,idx[i]]>0),which(connMat[idx[i],]>0)))
	if(length(coi) > 0){
		#create sub mat
		mati <- as(t(replicate(length(coi), matches[idx[i],])),"dgCMatrix")
		#add it to sub mat of connected peaks
		matches2[coi,,drop=FALSE] <- matches2[coi,,drop=FALSE] + mati
	}
}
diff <- Matrix::colSums(matches2) - Matrix::colSums(matches)
#print(diff)

#Make Annotation SE
message("INFO : Make Annotation SummarizedExperiment")
anno_ix <- SummarizedExperiment(assays=SimpleList(motifMatches=matches2), rowRanges=rowRanges(se))

#Compute Deviations
message("INFO : Compute Deviations")
assayNames(se)="counts"
dev <- computeDeviations(se, anno_ix)

#compute variability
message("INFO : compute variability")
metadata(dev)$Variability <- computeVariability(dev)

#add matches
metadata(dev)$gwasMatches <- anno_ix

#save output
message("INFO : Save Result")
saveRDS(dev, file.path(args$outdir,"GWAS-Co-Accessibility-chromVAR-Summarized-Experiment.rds"))
