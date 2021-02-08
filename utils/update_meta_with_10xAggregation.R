library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)
library(RColorBrewer)
library(viridis)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
		    default=NULL,
                    help="the path of project saved")

#parser$add_argument("--column",
#                    type="character",
#                    default="Sample",
#		    help="the column in cellcol to be exported")

parser$add_argument("--aggr_csv",
                    type="character",
                    default=NULL,
		    help="the aggregation.csv  from cellranger aggr or reanalyze")

#parser$add_argument("--outdir",
#                    type="character",
#                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()

options(stringsAsFactors=FALSE)
#if(!dir.exists(args$outdir)){
#        dir.create(args$outdir,recursive=TRUE)
#}



print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# Loading ArrowFiles")
projHeme<-loadArchRProject(args$project)


print("# Get cellcolData from ArchR Project")
cellCol=getCellColData(projHeme)
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
projHeme=addCellColData(projHeme,data=DATA$barcode,cells =DATA$cells, name = "Barcode",force=TRUE)
if("status"%in%str_to_lower(colnames(DATA))){
	projHeme=addCellColData(projHeme,data=DATA$status,cells =DATA$cells, name = "Status",force=TRUE)
}
DATA=DATA[rownames(cellCol),]
#cellCol=cbind(cellCol,DATA)
print("# Update ArchR object")
path=getOutputDirectory(projHeme)
saveRDS(projHeme,file.path(path,"Save-ArchR-Project.rds"))





