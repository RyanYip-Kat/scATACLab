library(Seurat)
library(future)
library(dplyr)

source("/home/ye/Work/R/scATAC/ArchR/plotter/ArchRHeatmap.R")

N=30
groupBy="RNA_snn_res.1"
Node="scale.data"
metaData=seurat@meta.data

plan("multiprocess", workers = 16)
markers=FindAllMarkers(object=seurat,
                       logfc.threshold = 0.25,
                       test.use = "wilcox",
                       slot="data")
topN=markers%>%group_by(cluster)%>%top_n(N,wt=avg_log2FC)
mat=GetAssayData(seurat,slot=Node)
genes=topN$gene[topN$gene%in%rownames(mat)]
rm.genes=topN$gene[!topN$gene%in%rownames(mat)]

if(length(rm.genes)>0){
	txt=paste(rm.genes,collapse=",")
	cat(sprintf("INFO : Genes :  %s\t\tnot in slot matrix\n",txt)) 
}


################  reorder cells
#groups=unique(Idents(seurat))
groups=levels(seurat)
cellSets=unlist(lapply(groups,function(x)return(WhichCells(seurat,idents=x))))
colData=metaData[cellSets,groupBy,drop=FALSE]
m=mat[genes,cellSets]


##############
if(Node!="scale.data"){
        limits=c(-2.5,2.5)
        scale-TRUE
}else{
	limits=c(min(m),max(m))
	scale=FALSE
}



p=ArchRHeatmap(m,colData=colData,
	       clusterRows=F,
	       clusterCols=F,
	       labelRows=T,
	       draw=FALSE,
	       limits=limits,
	       scale=scale)
