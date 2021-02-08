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
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--save",
                    action="store_true",
                    default=FALSE,
		    help="whether to save subset ArchR Project")

parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")


parser$add_argument("--tss_frag",
                    type="character",
                    default=NULL,
		    help="csv file with column : sample,tss,frag")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
print(paste0("### Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads)

########################
max_frags=50000
max_tss=15

#########################
print("### Loading Project")
projHeme<-loadArchRProject(args$project)
metric_df=as.data.frame(getCellColData(projHeme))
metric_df=subset(metric_df,nFrags < max_frags & TSSEnrichment < max_tss)

tss_frag_df=read.csv(args$tss_frag,stringsAsFactors=FALSE)
colnames(tss_frag_df)=c("Sample","Tss","Frag")
rownames(tss_frag_df)=tss_frag_df$Sample

sample_lists=unique(metric_df$Sample)
cell_lists=c()
for(i in 1:length(sample_lists)){
	sam=sample_lists[i]
	metric_df_subset=subset(metric_df,str_detect(Sample,sam))
	print("------------------")
	print(paste0("Pre QC : #nCell of : ",sam," is : ",nrow(metric_df_subset)))
	tss_flag=subset(tss_frag_df,str_detect(Sample,sam))
	tss=as.numeric(tss_flag$Tss)
	frag=as.numeric(tss_flag$Frag)
	print(paste0("Filter Sample :",sam," with Tss : ",tss," Frag : ",frag))
	metric_df_subset=subset(metric_df_subset,TSSEnrichment > tss & nFrags >frag)
	cells=rownames(metric_df_subset)
	cell_lists=c(cell_lists,cells)
	print(paste0("After QC : #nCell of : ",sam," is : ",length(cells)))
	print("------------------")
}

print(paste0("After QC : #nCell of total is  : ",length(cell_lists)))
write.table(cell_lists,file.path(args$outdir,"qc_cells.csv"),quote=FALSE,row.names=FALSE,sep=",",col.names=FALSE)
if(args$save){
	print("### Save subset ArchRProject")
	subsetArchRProject(projHeme,cells=cell_lists,outputDirectory=file.path(args$outdir,"ArchRSubset-with-sample"))
}



