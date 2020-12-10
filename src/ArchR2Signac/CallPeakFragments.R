library(Seurat)
library(GenomicRanges)
library(Signac)
library(argparse)
library(stringr)
library(future)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")

parser$add_argument("--assay",
                    type="character",
                    default="ATAC",
		    choices=c("peaks","ATAC","atac"),
                    help="assay used")

parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")

args <- parser$parse_args()

###############################
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}


###############################  Configure
Genome=args$genome

if(Genome=="hg38"){
        library(EnsDb.Hsapiens.v86)
        annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
}else{
        library(EnsDb.Mmusculus.v79)
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
}

# change to UCSC style since the data was mapped to hg19
cat(sprintf("INFO : Annotation with [ %s ]\n",Genome))
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- Genome

############################## Configure
plan("multiprocess", workers = 16)
outDir=dirname(args$seurat)
callpeakDir=file.path(outDir,"CallPeak")

makedir(callpeakDir)
############################## Loading dataset
message("INFO : Loading dataset")
seurat=readRDS(args$seurat)
DefaultAssay(seurat)=args$assay
fragments_list=Fragments(seurat)
fragments_name=names(fragments_list)
peak_list=list()
gr_list=GRangesList()
############################## callpeaks
for(i in seq_along(fragments_name)){
	name=fragments_name[i]
	outPath=file.path(callpeakDir,name)
	makedir(outPath)

	peakName=file.path(outPath,paste0(name,"_peaks.rds"))
	if(!file.exists(peakName)){
		cat(sprintf("INFO : callpeak with [ %d of %d ] --- [ %s ]\n",i,length(fragments_name),name))
		peaks <- CallPeaks(object =fragments_list[[name]] ,
			   macs2.path = findMacs2(),
			   outdir=outPath,
			   cleanup=FALSE,
			   effective.genome.size=2.7e+09)
		message("INFO : Save callpeak result ...")
	        saveRDS(peaks,file.path(outPath,paste0(name,"_peaks.rds")))

	}else{
		cat(sprintf("INFO : Loading peaks [ %d of %d ] --- [ %s ]\n",i,length(fragments_name),peakName))
		peaks=readRDS(peakName)
	}
	peak_list[[name]]=peaks
	gr_list[[name]]=peaks
}

###################################
combined.peaks=UnifyPeaks(gr_list)
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

message("INFO : get callpeak counts ...")
macs2_counts=FeatureMatrix(fragments=Fragments(seurat),features=combined.peaks,cells = colnames(seurat))

message("INFO : Create Object with callpeak counts ...")
metadata=seurat@meta.data
peaks.assay=CreateChromatinAssay(counts=macs2_counts,
                                  fragments=Fragments(seurat),
                                  genome =Genome,
                                  sep = c("-", "-"),
                                  annotation=annotations,
                                  validate.fragments = FALSE)
seurat.peaks <- CreateSeuratObject(peaks.assay, assay = "peaks",meta.data = metadata)
Annotation(seurat.peaks) <- annotations
seurat.peaks[["ATAC"]]=seurat[["ATAC"]]
message("INFO : Save Object ...")
saveRDS(seurat.peaks,file.path(outDir,"seurat-peaks.rds"))
message("INFO : Done!")

