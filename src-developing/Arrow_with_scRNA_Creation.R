library(argparse)
library(stringr)
library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
library(SummarizedExperiment)
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--files",
                    type="character",
                    default=NULL,
		    help="the csv file,include(arrow files,scRNA *.h5 files,Status,..")


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



args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print(paste0("Setting default genome  with :",args$genome))
addArchRGenome(args$genome)

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads)

inputFiles<-read.table(args$files,sep=",",stringsAsFactors=FALSE)
Samples=as.character(inputFiles$V1)  # sample name
ArrowFiles=as.character(inputFiles$V2)     # *.arrow
H5CountsFiles=as.character(inputFiles$V3)  # filtered_feature_bc_matrix.h5 
Seurats=as.character(inputFiles$V4)
print(head(inputFiles))
stopifnot(length(ArrowFiles)==length(H5CountsFiles))

print("### Loading Seurats")
seurat_cells=list()
for(i in 1:length(Seurats)){
	seurat=readRDS(Seurats[i])
	#counts=GetAssayData(seurat,"counts")
	cells=colnames(seurat)
	cells=unlist(lapply(cells,function(cell){
				    cell=str_split(cell,"-")[[1]][1]
				    return(paste0(cell,"-1"))})
	)
	cells=paste(Samples[i],cells,sep="#")
	print(head(cells))
	seurat_cells[[i]]=cells
}
print("### Import scRNA")
seRNA <- import10xFeatureMatrix(
    input = H5CountsFiles,
    names = Samples
)

seRNA <- Reduce("cbind", seRNA)
#print("### Subset seRNA from Seurat")
#for(i in 1:length(seurat_cells)){
#	cells=seurat_cells[[i]]
#	se=seRNA[[i]]
#	seRNA[[i]]=se[,cells]
#}

#counts_list=list()
#for(i in seq_along(H5CountsFiles)){
#	print(paste0("### ","Sample :",Samples[i],"  Loading gene matrix from :",H5CountsFiles[i]))
#	count=Read10X_h5(H5CountsFiles[i])
#	Names=unlist(lapply(colnames(count),function(x)return(str_split(x,"-")[[1]][1])))
#	Names=paste(Names,"1",sep="-")
#	cellNames=paste(Samples[i],Names,sep="#")
#	colnames(count)=cellNames
#	counts_list[[i]]=count
#}
#counts=do.call(cbind,counts_list)
#features=rownames(counts)
#se <- SummarizedExperiment(assays = SimpleList(counts = counts), rowData = features)


print("### Create ArchRProject")
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
print(head(getCellNames(proj)))
print("### addGeneExpressionMatrix from seRNA")
proj=addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

print("### addIterativeLSI")
proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(
      resolution = 0.2,
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    name = "LSI_ATAC"
)

print("## Combined ATAC and RNA Dims")
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")


#print("# Run harmony ")
#proj <- addHarmony(
#    ArchRProj = proj,
#    reducedDims ="LSI_Combined",
#    name = "Harmony",
#    groupBy = "Sample",
#    force = TRUE
#)

print("### add ATAC UMAPs")
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
print("### add RNA UMAPs")
proj <- addUMAP(proj, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)

print("### add combined UMAP")
proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)

print("### Add Clusters")
proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.8, force = TRUE)

print("### Save Project")
saveArchRProject(ArchRProj = proj, outputDirectory =file.path(args$outdir,"Save-ProjHeme-seRNA-scATAC"), load = TRUE)
