library(ArchR)
library(stringr)
library(dplyr)
library(argparse)


projHeme=loadArchRProject("/home/ye/Work/R/scATAC/ArchR/ArchR_Result/CEpi_Mouse/Save-ProjHeme-P2G/")
peakSet=getPeakSet(projHeme)

N=20
p2g <- metadata(projHeme@peakSet)$Peak2GeneLinks
P2GMatrix=plotPeak2GeneHeatmap(projHeme,groupBy="Labels",returnMatrices=TRUE,k=N)

#mATAC <- readRDS(metadata(p2g)$seATAC)
#mRNA <- readRDS(metadata(p2g)$seRNA)
#p2g$peak <- paste0(rowRanges(mATAC))
#p2g$gene <- rowData(mRNA)$name


############################
Peak2GeneLinks=P2GMatrix[["Peak2GeneLinks"]]
ATAC=P2GMatrix[["ATAC"]]
RNA=P2GMatrix[["RNA"]]

kmeansId=ATAC[["kmeansId"]]
DF=data.frame("peak"=Peak2GeneLinks$peak,"gene"=Peak2GeneLinks$gene,"kid"=as.factor(kmeansId))
rownames(DF)=rownames(Peak2GeneLinks)


DF_list=list()
for(id in unique(DF$kid)){
	cat(sprintf("INFO : select [ %d ] gene in [ %s ] kmeansId\n",N,id))
	df=subset(DF,kid==id)
	n=N
	if(N>nrow(df)){
		   n=nrow(df)
	}
	x=df[1:n,,drop=FALSE]
	x$rowName=rownames(x)
	DF_list[[id]]=x
}

topDF=do.call(rbind,DF_list)
rownames(topDF)=topDF$rowName

topDF$chr=unlist(lapply(topDF$peak,function(x)str_split(x,":")[[1]][1]))
topDF$start=unlist(lapply(topDF$peak,function(x)str_split(x,":|-")[[1]][2]))
topDF$end=unlist(lapply(topDF$peak,function(x)str_split(x,":|-")[[1]][3]))
topDFGR=makeGRangesFromDataFrame(topDF,keep.extra.columns=T)


###########################
idx=findOverlaps(topDFGR,peakSet)
topDFGR=topDFGR[idx@from,]
mcols(topDFGR)$nearstGene=mcols(peakSet[idx@to,])$nearestGene
saveRDS(topDFGR,"topDFGR.rds")
