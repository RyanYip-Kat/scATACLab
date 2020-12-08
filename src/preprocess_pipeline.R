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
                    default="")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--genome",
                    type="character",
                    default="hg38",
		    choices=c("hg38","hg19","mm9","mm10"),
		    help="which genome to be used")

parser$add_argument("--num_threads",
                    type="integer",
                    default=4,
                    help="number of threads for command")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

inputFiles<-read.table(args$inputFiles,sep=",",stringsAsFactors=FALSE)
files=as.character(inputFiles$V1)
filenames=as.character(inputFiles$V2)
barcode_csv=unlist(lapply(files,function(file){return(file.path(dirname(file),"analysis/tsne/2_components/projection.csv"))}))

print("# Get valid barcode")
csvFiles=unlist(lapply(1:length(barcode_csv),function(i){
               csv=barcode_csv[i]
               d=read.csv(csv,stringsAsFactors=FALSE,header=TRUE)
               d=subset(d,select=Barcode)
               colnames(d)="barcode"
               d$cell_id=1:nrow(d)
               filename=file.path(args$outdir,paste0(filenames[i],"_validBarcode.csv"))
               write.table(d,filename,sep=",",quote=F,row.names=F)
               return(filename)
}))
sampleName=unlist(lapply(csvFiles,function(file){return(str_split(basename(file),"_")[[1]][1])}))

print("----------")
print(files)
print("##########")
print(sampleName)
print("##########")
print(csvFiles)
print("----------")

print("# Get valid barcode from cellranger counts")
barcodeList=getValidBarcodes(csvFiles,sampleName)
print(paste0("Setting default genome  with :",args$genome))
addArchRGenome(args$genome)

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads)

print("# Create ArrowFiles")
ArrowFiles <- createArrowFiles(
  inputFiles = files,
  sampleNames =filenames,
  validBarcodes=barcodeList,  # with the barcode from cellranger
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
)



doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1,
)

projHeme <- ArchRProject(
  ArrowFiles = ArrowFiles,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

print("# Plot : Log10 Unique Fragments vs TSS Enrichment")
df <- getCellColData(projHeme, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")


plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme, addDOC = FALSE)

p1 <- plotGroups(
    ArchRProj = projHeme,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
   )

p2 <- plotGroups(
    ArchRProj = projHeme,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

p3 <- plotGroups(
    ArchRProj = projHeme,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "ridges"
   )

p4 <- plotGroups(
    ArchRProj = projHeme,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projHeme, addDOC = FALSE, width = 12, height = 12)

p1 <- plotFragmentSizes(ArchRProj = projHeme)
p2 <- plotTSSEnrichment(ArchRProj = projHeme)

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projHeme, addDOC = FALSE, width = 8, height = 8)

saveArchRProject(ArchRProj = projHeme, outputDirectory =file.path(args$outdir,"Save-ProjHeme"), load = TRUE)
