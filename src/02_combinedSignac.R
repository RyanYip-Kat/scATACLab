library(stringr)
library(DropletUtils)
library(Seurat)
library(Signac)
library(argparse)
library(GenomicRanges)
library(future)

source("/home/ye/Work/R/scATAC/Signac/src/helper_functions.R")
parser <- ArgumentParser(description='Process to combined mutiple Signac Objects')
parser$add_argument("--meta",
                    type="character",
                    default=NULL,
                    help="meta *.csv file(sample,*.rds,...)")


parser$add_argument("--outdir",
                    type="character",
                    default="./Result")

parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")


args <- parser$parse_args()
############################### Configure
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

plan("multiprocess", workers = 16)
options(future.globals.maxSize = 960000 * 1024^2) # for 50 Gb RAM
outDir=args$outdir
Genome=args$genome
makedir(outDir)

if(Genome=="hg38"){
        library(EnsDb.Hsapiens.v86)
        annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
	annoFile="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
}else{
        library(EnsDb.Mmusculus.v79)
        annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
	annoFile="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-mm10-3.0.0/genes/genes.gtf"
}
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- Genome


############################### Loading datasets
DF=read.csv(args$meta,sep=",",stringsAsFactors=FALSE,header=FALSE)
Names=as.character(DF$V1)
Files=as.character(DF$V2)

###############################
message("INFO : Loading dataset")
sep="_"
fragment_list=list()
seurat_list=list()
cell_list=list()
for(i in seq_along(Files)){
	    name=Names[i]
	    file=Files[i]
	    cat(sprintf("INFO : [ time : %s ] --- [ %d of %d] --- [ %s from %s] \n",Sys.time(),i,length(Files),name,file))
	    seurat=readRDS(file)
	    counts=GetAssayData(seurat,"counts")
	    colnames(counts)=paste(name,colnames(counts),sep=sep)
	    fragment=Fragments(seurat)[[1]]@path
	    fragment=ReadFragments(fragment,name,sep)
	    seurat_assay <- CreateChromatinAssay(
                                             counts = counts,
                                             sep = c("-", "-"),
                                             genome =Genome,
                                             fragments = fragment)
	    seurat <- CreateSeuratObject(
                                     counts = seurat_assay,
                                     assay = 'ATAC',
                                     project = 'scATAC')

	    seurat_list[[name]]=seurat
	    fragment_list[[name]]=fragment
	    cell_list[[name]]=Cells(seurat)
}

saveRDS(cell_list,file.path(outDir,"cell_list.rds"))
###############################
combined.peaks=UnifyPeaks(seurat_list)

####################################   Merge all Object
combined_list=list()
fragment_assay_list=list()
meta_list=list()
message("INFO : Make Feature Matrix ...")
i=1
for(name in names(fragment_list)){
        cat(sprintf("INFO : [ time : %s ] --- [ %d of %d ]\n",Sys.time(),i,length(fragment_list)))
        cells=cell_list[[name]]
        fragment=fragment_list[[name]]

        message("INFO : Create FragmentObject")
        frags=CreateFragmentObject(path=fragment,cells=cells,validate.fragments=FALSE)
        message("INFO : Make FeatureMatrix ")
        Fragcount<- FeatureMatrix(
                                  fragments = frags,
                                  features = combined.peaks,
                                  cells =cells)

        message("INFO : Create Seurat Object ...")
        Fassay <- CreateChromatinAssay(Fragcount, fragments = frags)
        seurat <- CreateSeuratObject(Fassay, assay = "ATAC")
	fragment_assay_list[[name]]=Fragments(seurat)[[1]]
        #seurat=AddMetaData(seurat,meta_list[[name]])
        seurat$dataset=name
	seurat$Cells=Cells(seurat)
        cat(sprintf("INFO : [ time : %s ] --- Object Size  : [ %d , %d ] \n",Sys.time(),nrow(seurat),ncol(seurat)))
        combined_list[[name]]=seurat
	meta_list[[name]]=seurat@meta.data
        i=i+1
}

message("INFO : Merge Object ...")
###################################
message("INFO : Merge Count lists ...")
#seurat=Reduce(merge,seurat_list)
count_list=lapply(combined_list,function(seurat){return(GetAssayData(seurat,"counts",assay="ATAC"))})
counts=do.call(cbind,count_list)
##################################### Create Combined Object
metadata=do.call(rbind,meta_list)

#seurat=merge(x=combined_list[[1]],
#             y=combined_list[-1])
             #add.cell.ids=names(combined_list))  # add name+"_" prefix
#seurat=Reduce(merge,combined_list)
##################################### Create Combined Object
message("INFO : Make Combined Object ...")
#counts=GetAssayData(seurat,"counts",assay="ATAC")
fragments=fragment_assay_list
#metadata=seurat@meta.data

peaks.assay=CreateChromatinAssay(counts=counts,
				  fragments=fragments,
				  genome =Genome,
				  sep = c("-", "-"),
				  annotation=annotations,
				  validate.fragments = FALSE)
seurat <- CreateSeuratObject(peaks.assay, assay = "ATAC")
rownames(metadata)=metadata$Cells
metadata=metadata[Cells(seurat),]
seurat=AddMetaData(seurat,metadata)
Annotation(seurat) <- annotations


##################################### bug not fixed
#message("INFO : Add RNA Slot with  Gene Activity ...")
#message("INFO : Get GeneActivity")
#gene.activities <- GeneActivity(seurat)
#activity.matrix <- CreateGeneActivityMatrix(peak.matrix = counts, annotation.file = annoFile,
#    seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE,keep.sparse=TRUE)

#seurat[['RNA']] <- CreateAssayObject(counts = activity.matrix)
#seurat <- NormalizeData(
#  object = seurat,
#  assay = 'RNA',
#  normalization.method = 'LogNormalize',
#  scale.factor = median(object$nCount_RNA)
#)

message("INFO : Save Combined Object")
saveRDS(seurat,file.path(outDir,"combined_seurat.rds"))
message("INFO : Done!")

