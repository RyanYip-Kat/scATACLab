library(ArchR)
library(stringr)
library(DropletUtils)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(argparse)
library(GenomicRanges)
library(future)


source("/home/ye/Work/BioAligment/SNP/Shi/scripts/help_functions.R")
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="the project path  of ArchR")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--column",
                    type="character",
                    default=NULL,
                    help="celltype column's name")

parser$add_argument("--subset",
                    type="character",
                    nargs="+",
                    default=NULL,
                    help="subset's for subsetting")

parser$add_argument("--genome",
                    type="character",
                    default="hg38",
                    choices=c("hg38","mm10"),
                    help="which genome to be used")

args <- parser$parse_args()

#############################
plan("multiprocess", workers = 24)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM

outDir=args$outdir
fragmentDir=file.path(outDir,"fragments")
bedDir=file.path(outDir,"bedFile")
makedir(outDir)
makedir(fragmentDir)
makedir(bedDir)

project=args$project

bgzip="/home/ye/anaconda3/envs/BulkBio/bin/bgzip"
tabix="/home/ye/anaconda3/envs/BulkBio/bin/tabix"

#############################
ZipAndBuildTbi=function(gzipFile){
	fileName=file.path(dirname(gzipFile),paste0("sorted-",basename(gzipFile)))
	cmd=paste0("sort -k1,1d -k2,2n -k3,3n"," ",gzipFile," > ",fileName)
	system(cmd)
	cmd=paste0(bgzip," ",fileName)
	message("INFO : Zip File ...")
	system(cmd)

	cmd=paste0(tabix," ","--preset=bed"," ",fileName,".gz")
	message("INFO : Build tabix Index  ...")
	system(cmd)
	system(paste0("rm ",gzipFile))
	return(paste0(fileName,".gz"))
}


message("INFO : Loading dataset")
addArchRThreads(threads=16)
projHeme=loadArchRProject(project)
metadata=as.data.frame(getCellColData(projHeme))
cells=getCellNames(projHeme)

projHeme=loadArchRProject(args$project)
metadata=as.data.frame(getCellColData(projHeme))
cells=getCellNames(projHeme)

############################## subset

if(!is.null(args$column) & !is.null(args$subset)){
	assertthat::assert_that(args$column%in%colnames(metadata))
	cat(sprintf("INFO : Subset from [ %s ] with : [ %s ]",args$column,paste(args$subset),collapse=","))
        target=metadata[[args$column]]

        stopifnot(args$subset%in%unique(target))
        idxPass= which(target%in%args$subset)
        cellPass=cells[idxPass]
        projHeme=subsetCells(projHeme,cellNames=cellPass)
	metadata=as.data.frame(getCellColData(projHeme))
        cells=getCellNames(projHeme)
}



############################### Whether to save scATAC peaks matrix
peakFile=file.path(outDir,"scATAC-PeakMatrix.rds")
if(!file.exists(peakFile)){
	############################### get Peakset and PeakMarix
        message("INFO : getPeakSets")
        peaks=getPeakSet(projHeme)
        #cells=unlist(lapply(cells,function(cell)return(str_split(cell,"#")[[1]][2])))

        message("INFO : Get PeakMatrix")
        se=getMatrixFromProject(projHeme,useMatrix="PeakMatrix")
        rowData(se)=peaks
        all.out=cbind(data.frame(chr=seqnames(peaks)),as.data.frame(ranges(peaks)))
        peak=paste0(all.out$chr,":",all.out$start,"-",all.out$end)
        rownames(se)=peak
        #colnames(se)=cells
        #counts=assay(se)

	message("INFO : Save PeakMatrix")
	saveRDS(se,peakFile)

}else{
	message("INFO : loading PeakMatrix ...")
	se=readRDS(peakFile)
}

############################### Whether to save scATAC GeneScore matrix
GAFile=file.path(outDir,"scATAC-GA-Matrix.rds")
if(!file.exists(GAFile)){
	#######################   get Gene Score Matrix
        message("INFO : Get Score GeneMatrix")
        seRNA=getMatrixFromProject(projHeme,useMatrix = "GeneScoreMatrix")
        features=getFeatures(projHeme,useMatrix = "GeneScoreMatrix")
        rownames(seRNA)=features

	message("INFO : Save Gene Score Matrix")
	saveRDS(seRNA,GAFile)
}else{
        message("INFO : loading Gene Score Matrix ...")
        seRNA=readRDS(GAFile)
}


############################### get Fragment Files
metadata=as.data.frame(getCellColData(projHeme))
cellNames=getCellNames(projHeme)
metadata$Cells=cellNames

arrowFiles=getArrowFiles(projHeme)
Donors=names(arrowFiles)
keepChrs=paste("chr",c(1:22,"X","Y"),sep="")

############################### list objects
fragment_list=list()
count_list=list()
score_list=list()
meta_list=list()
bed_list=list()
cell_list=list()

###############################
message("INFO : get Fragment Files and counts matrix ...")
counts=assay(se)
score=assay(seRNA)
for(i in seq_along(arrowFiles)){
	donor=Donors[i]
	cat(sprintf("INFO : [ %s of %s ] : [ %s ] --- [ %s ]\n",i,length(arrowFiles),donor,arrowFiles[i]))
	keepCells=metadata$Cells[metadata$Sample==donor]  # sample#barcode

	zipFrag=file.path(fragmentDir,paste0("sorted-",donor,"_fragment.tsv.gz"))
	if(!file.exists(zipFrag)){
		fragment=getFragmentsFromArrow(ArrowFile=arrowFiles[i],chr=keepChrs,cellNames=keepCells)
	        fragment=as.data.frame(fragment)
	        #cell=unlist(lapply(fragment$RG,function(cell)return(str_split(cell,"#")[[1]][2])))
	        #fragment$cells=cell

		outFrag=file.path(fragmentDir,paste0(donor,"_fragment.tsv"))
		cat(sprintf("INFO : Write fragment into : [ %s ] \n",outFrag)) 
	        write.table(fragment[,c("seqnames","start","end","RG")],outFrag,sep="\t",row.names=F,quote=F,col.names=F)

	        zipFrag=ZipAndBuildTbi(outFrag)
	        fragment_list[[donor]]=zipFrag
	}else{
		fragment_list[[donor]]=zipFrag
	}

	count=counts[,keepCells]
	GA=score[,keepCells]
	meta=metadata[keepCells,]
	#cell=unlist(lapply(keepCells,function(cell)return(str_split(cell,"#")[[1]][2])))
	cell=keepCells
	colnames(count)=cell
	colnames(GA)=cell
	rownames(meta)=cell
	cat(sprintf("INFO : [ %s ] Matrix's Size is [ %d,%d ]\n",donor,dim(count)[1],dim(count)[2]))
	count_list[[donor]]=count
	score_list[[donor]]=GA
	meta_list[[donor]]=meta
	cell_list[[donor]]=cell
}

saveRDS(fragment_list,file.path(outDir,"fragment_list.rds"))
#saveRDS(meta_list,file.path(outDir,"meta_list.rds"))
saveRDS(score,file.path(outDir,"gA.rds"))
###############################  Create Signac Object
seurat_list=list()
rna_list=list()
message("INFO : Create Signac Object ...")
for(name in names(count_list)){
	cat(sprintf("INFO : Create Signac Object : %s \n",name))
	count=count_list[[name]]
	fragment=fragment_list[[name]]
	CAssay=CreateChromatinAssay(counts=count,
				   sep = c(":", "-"),
				   genome="hg38",
				   fragments=fragment,validate.fragments=F)
	seurat=CreateSeuratObject(CAssay,assay = 'peaks',project = 'ATAC')
	seurat=AddMetaData(seurat,meta_list[[name]])
	seurat$dataset=name
	cat(sprintf("INFO :  Object Size  : [ %d , %d ] \n",nrow(seurat),ncol(seurat)))

	cat(sprintf("INFO : Add : %s Gene Score \n",name))
	archrGA <- CreateAssayObject(counts =score_list[[name]] )
	#seurat[['archrGA']] <- archrGA
        rna_list[[name]]=archrGA	
	#seurat <- NormalizeData(
	#		      	object = seurat,
        #                        assay = 'archrGA',
        #                        normalization.method = 'LogNormalize')
	seurat_list[[name]]=seurat
}

saveRDS(rna_list,file.path(outDir,"GA_list.rds"))

###############################  Annotation Peaks 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- args$genome #"hg38"

message("INFO : Annotate Signac Object")
bed_list=list()
cell_list=list()
for(name in names(seurat_list)){
	#cat(sprintf("INFO : Annotation  [ %d of %d ] \n",i,length(seurat_list)))
	seurat=seurat_list[[name]]
	Annotation(seurat) <- annotations
	seurat_list[[name]]=seurat
	cell_list[[name]]=Cells(seurat)

	peakName=file.path(bedDir,paste0(name,".bed"))
	if(!file.exists(peakName)){
		peaks=rownames(seurat)
	        Lists=lapply(peaks,function(peak){return(str_split(peak,":|-")[[1]])})
		chrom=unlist(lapply(Lists,function(x)x[1]))
		st=unlist(lapply(Lists,function(x)x[2]))
		ed=unlist(lapply(Lists,function(x)x[3]))
		peakDF=data.frame("chr"=chrom,"start"=st,"end"=ed)
		write.table(peakDF,peakName,sep="\t",row.names=F,col.names=F,quote=F)
	}
	bed_list[[name]]=peakName

}

message("INFO : Save ..")
saveRDS(bed_list,file.path(outDir,"bed_list.rds"))
saveRDS(cell_list,file.path(outDir,"cell_list.rds"))
############################### Combined Peaks
makePeakGR<-function(bedFile){
	peak= read.table(
			 file = bedFile,
			 col.names = c("chr", "start", "end"))
	gr=makeGRangesFromDataFrame(peak)
	return(gr)
}


GR_list=list()
for(name in names(bed_list)){
	gr=makePeakGR(bed_list[[name]])
	GR_list[[name]]=gr
}
saveRDS(GR_list,file.path(outDir,"gr-list.rds"))


###################################
combined.peaks=UnifyPeaks(seurat_list)
####################################   Merge all Object
combined_list=list()
fragment_assay_list=list()
meta_list=list()
message("INFO : Make Feature Matrix ...")
i=1
for(name in names(fragment_list)){
	cat(sprintf("INFO : [ %d of %d ]\n",i,length(fragment_list)))
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
	seurat=AddMetaData(seurat,meta_list[[name]])
        seurat$dataset=name
        cat(sprintf("INFO :  Object Size  : [ %d , %d ] \n",nrow(seurat),ncol(seurat)))
	#seurat[["archrGA"]]=rna_list[[name]]
	#seurat <- NormalizeData(
        #                        object = seurat,
        #                        assay = 'archrGA',
        #                        normalization.method = 'LogNormalize')
	combined_list[[name]]=seurat
	meta_list[[name]]=seurat@meta.data
	fragment_assay_list[[name]]=Fragments(seurat)[[1]]
	i=i+1
}

saveRDS(meta_list,file.path(outDir,"meta_list.rds"))
saveRDS(combined_list,file.path(outDir,"seurat_list.rds"))
saveRDS(fragment_assay_list,file.path(outDir,"fragment_assay_list.rds"))


############################## get reduced slots
message("INFO : get ArchR ReduceDims for Signac Object")
ReducedDimNames=getReducedDimsNames(projHeme)
EmbeddingNames=getEmbeddingsNames(projHeme)
EmbeddingNames=setdiff(EmbeddingNames,"UMAPHarmony") # not include "UMAPHarmony"
#cellNames=str_replace(cellNames,"#","_")

ReducedDimDir=file.path(outDir,"ReducedDims")
EmbeddingDir=file.path(outDir,"Embeddings")
makedir(ReducedDimDir)
makedir(EmbeddingDir)
for(name in EmbeddingNames){
	cat(sprintf("INFO : Extract [ %s ]\n",name))
        mat=as.matrix(getEmbedding(ArchRProj =projHeme, embedding =name, returnDF = TRUE))
	#rownames(mat)=cellNames
	fileName=file.path(EmbeddingDir,paste0(name,".rds"))
	saveRDS(mat,fileName)
}

for(name in ReducedDimNames ){
	cat(sprintf("INFO : Extract [ %s ]\n",name))
        mat=getReducedDims(projHeme,reducedDims=name,returnMatrix=T)
	#rownames(mat)=cellNames
	fileName=file.path(ReducedDimDir,paste0(name,".rds"))
	saveRDS(mat,fileName)
}




