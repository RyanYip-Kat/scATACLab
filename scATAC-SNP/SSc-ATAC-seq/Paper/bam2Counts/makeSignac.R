library(Seurat)
library(GenomicRanges)
library(Signac)
library(argparse)
library(stringr)
library(future)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--meta",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")

parser$add_argument("--assay",
                    type="character",
                    default="ATAC",
                    choices=c("peaks","ATAC","atac"),
                    help="assay used")

parser$add_argument("--unionPeak",
		    type="character",
		    default=NULL,
		    help="union peak *.rds file")

parser$add_argument("--outdir",
                    type="character",
                    default="Results")

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
outDir=args$outdir

###############################
DF=read.csv(args$meta,sep=",",header=FALSE)
fragment_list=as.character(DF$V2)
name_list=as.character(DF$V1)
combined.peaks=readRDS(args$unionPeak)

seurat_list=list()
for(i in seq_along(fragment_list)){
	name=name_list[i]
	cat(sprintf("INFO : Create Fragment Object [ %d of %d ] --- [ %s ]\n",i,length(fragment_list),name))
	frag=CreateFragmentObject(path=fragment_list[i],tolerance = 0.5)
        counts <- FeatureMatrix(
				fragments = frag,
				features = combined.peaks)
	fassay <- CreateChromatinAssay(counts, 
				       fragments = frag,
				       genome =Genome,
				       annotation=annotations,
				       sep = c("-", "-"),
				       validate.fragments = FALSE)
	seurat <- CreateSeuratObject(fassay, assay = args$assay,min.cells=10,min.features=200)
	seurat$dataset=name
	cat(sprintf("INFO : [ %s ] Seurat Object Size -- [ %d , %d ]\n",name,nrow(seurat),ncol(seurat)))
	seurat_list[[name]]=seurat
}

message("INFO : Save Object")
saveRDS(seurat_list,file.path(outDir,"seurat_list.rds"))
message("INFO : Done!")
