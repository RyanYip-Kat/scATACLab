library(argparse)
library(stringr)
library(jsonlite)
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

parser$add_argument("--aggr_csv",
                    type="character",
                    default=NULL,
                    help="the aggregation.csv  from cellranger aggr or reanalyze")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


inputFiles<-read.table(args$inputFiles,sep=",",stringsAsFactors=FALSE)
files=as.character(inputFiles$V1)
filenames=as.character(inputFiles$V2)
barcode_json=as.character(inputFiles$V3)# json file from clean_barcode_multiplets_1.1.py  and filter_barcode.py
print("# Get valid barcode")
csvFiles=unlist(lapply(1:length(barcode_json),function(i){
               json_file=barcode_json[i]
               d=as.data.frame(fromJSON(json_file))
               colnames(d)="barcode"
               d$cell_id=1:nrow(d)
	       filename=file.path(args$outdir,paste0(filenames[i],"_validBarcode.csv"))
               write.table(d,filename,sep=",",quote=F,row.names=F)
	       return(filename)
}))
sampleName=unlist(lapply(csvFiles,function(file){return(str_split(basename(file),"_")[[1]][1])}))
#barcodeList=list()
#for(i in 1:length(barcode_csv)){
#	csv=barcode_csv[i]
#	d=read.csv(csv,stringsAsFactors=FALSE,header=TRUE)
#	d=d[d$cell_id!="None",]
#	barcode=d$barcode
#	barcodeList[[filenames[i]]]=barcode
#}

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
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
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

if(args$rm_doublet){
	projHeme2 <- filterDoublets(projHeme)
}else{
	projHeme2=projHeme
}

print("# Get cellcolData from ArchR Project")
cellCol=getCellColData(projHeme2)
cellCol=as.data.frame(cellCol)

Sample=cellCol$Sample
cells=rownames(cellCol)

df=read.csv(args$aggr_csv,stringsAsFactors=F)
df$ident=1:nrow(df)
library_id=as.character(df$library_id)
ident=df$ident
barcode_list=list()
for(i in 1:length(library_id)){
        id=library_id[i]
        print(paste0("## Get barcode from : ",id))
        idx=str_detect(cells,id)
        cell_s=cells[idx]
        barcode=unlist(lapply(cell_s,function(s){return(str_split(s,"#")[[1]][2])}))
        barcode=unlist(lapply(barcode,function(s){
                                              x=str_split(s,"-")[[1]][1]
                                              x=paste0(x,"-",ident[i])
                                              return(x)}))

        data=data.frame("barcode"=barcode,row.names=cell_s)
        if("status"%in%colnames(df)){
                        data$status=df$status[i]
        }

        barcode_list[[i]]=data
}
DATA=do.call(rbind,barcode_list)
DATA$cells=rownames(DATA)
DATA=DATA[rownames(cellCol),]
projHeme2=addCellColData(projHeme2,data=as.character(DATA$barcode),cells =DATA$cells, name = "Barcode")
if("status"%in%str_to_lower(colnames(DATA))){
        projHeme2=addCellColData(projHeme2,data=as.character(DATA$status),cells =DATA$cells, name = "Status")
}

#cellCol=cbind(cellCol,DATA)
print("# Update ArchR object")

saveArchRProject(ArchRProj = projHeme2, outputDirectory =file.path(args$outdir,"Save-ProjHeme"), load = TRUE)
