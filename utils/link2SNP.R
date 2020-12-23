library(SummarizedExperiment)
library(stringr)
library(argparse)
library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)

featureName <- function(gr){
  paste(seqnames(gr),start(gr),end(gr),sep="_")
}



DF=read.table("tf2target/degs-hcbd-cdc-logfc0.01.csv",sep=",",stringsAsFactors=F,header=TRUE,row.names=1)   # from Seurat FindAllMarkers function
#DF=read.csv("tf2target/_markers.csv",sep=",",stringsAsFactors=F,header=TRUE)
Cols=unique(DF$cluster)
LFC=list()
Pval=list()
FDR=list()

#intersectGene_list=list()
#for(col in Cols){
#	df=subset(DF,cluster==col)
#	intersectGene_list[[col]]=df$names
#}

#intersectGenes=Reduce(intersect,intersectGene_list)
for(col in Cols){
	#df=subset(DF,cluster==col & names%in%intersectGenes)
	df=subset(DF,cluster==col)
	fc=df[,"avg_logFC",drop=FALSE]  #seurat
	rownames(fc)=df$gene
	#fc=df[,"logfoldchanges",drop=FALSE]
	#rownames(fc)=df$names
	colnames(fc)=col
	LFC[[col]]=fc

	pval=df[,"p_val",drop=FALSE]
	rownames(pval)=df$gene
	#pval=df[,"pvals",drop=FALSE]
	#rownames(pval)=df$names
	colnames(pval)=col
	Pval[[col]]=pval

	fdr=df[,"p_val_adj",drop=FALSE]
	rownames(fdr)=df$gene
	#fdr=df[,"pvals_adj",drop=FALSE]
	#rownames(fdr)=df$names
	colnames(fdr)=col
	FDR[[col]]=fdr
}

LFC_DF=do.call(cbind,LFC)
Pval_DF=do.call(cbind,Pval)
FDR_DF=do.call(cbind,FDR)

diffObj=SummarizedExperiment(assays=SimpleList("LFC"=LFC_DF,"FDR"=FDR_DF,"Pval"=Pval_DF))
colnames(diffObj)=Cols

diffATAC=readRDS("tf2target/degs-markerTest-hcbd-cdc.rds")
p2gLinks=readRDS("tf2target/Save_HCBD_P2G_Links.rds")$linksSig
matches=readRDS("tf2target/Matches.rds")

################################

rownames(matches) <- paste(seqnames(matches),start(matches),end(matches), sep = "_")
rownames(diffATAC) <- paste(seqnames(diffATAC),start(diffATAC),end(diffATAC), sep = "_")

#Make P2G Mats
names(p2gLinks) <- paste0("l",seq_along(p2gLinks))
dAp2g <- diffATAC[featureName(p2gLinks),colnames(diffObj)]
dRp2g <- diffObj[mcols(p2gLinks)$gene_name[mcols(p2gLinks)$gene_name%in%rownames(diffObj)],]
rownames(dAp2g) <- names(p2gLinks)
rownames(dRp2g) <- names(p2gLinks)[mcols(p2gLinks)$gene_name%in%rownames(diffObj)]
dAp2g=dAp2g[rownames(dRp2g),]

pdR=paste(seqnames(dAp2g),start(dAp2g),end(dAp2g), sep = "_")
matches=matches[pdR,]

#Identify Significant Peaks and Genes per MPAL Subpopulation
sigMatP2G <- abs(assays(dAp2g)[["Log2FC"]]) >= 0.1 & assays(dAp2g)[["Pval"]] <= 0.05 & abs(assays(dRp2g)[["LFC"]]) >= 0.1 & assays(dRp2g)[["Pval"]] <= 0.05

#Which have at least 1 MPAL subpopulation with a linked (within p2glinks) diff peak and diff gene
sigP2G <- p2gLinks[which(rowSums(sigMatP2G) > 0)]

#List of TFs to identify Targerts
tfs <-  colnames(matches)
#Identify Targets
t2gDF <- lapply(seq_along(tfs), function(x){
			tryCatch({
				#Subset Linked Peaks that contain TF Motif
				peaksx <- names(which(assay(matches[,tfs[x]])[,1]))
				#Subset sigMat by linked peaks with TF Motif
				sigMatP2Gx <- sigMatP2G[rownames(sigMatP2G) %in% peaksx,]
				#Subset links by linked peaks with TF Motif
				linksx <- sigP2G[featureName(sigP2G) %in% peaksx]
				#Figure out which samples are true for each link
                                sigMatP2Gx <- sigMatP2G[names(linksx),]

                                #Compute sum of MPAL subpopulations that are diff peak and diff peak for every peak connected to the gene
                                sigDFListx <- split(data.frame(sigMatP2Gx), mcols(linksx)$gene_name) %>% lapply(., colSums)

                                #Convert to Data Frame
                                sigDFx <- data.frame(Reduce("rbind",sigDFListx))
                                rownames(sigDFx) <- names(sigDFListx)

                                #Set to maximum of 1 for each MPAL subpop for all diff peak gene combo
                                sigDFx[sigDFx > 1] <- 1

                                #Compute number of positive MPAL subpopulations
                                nSubpop <- rowSums(sigDFx)
			        cat(sprintf("INFO : [ N = %d , P = %f ]\n",nSubpop,nSubpop / ncol(sigMatP2Gx)))      
                                #Summed R2 linkage for each linked positive diff peak gene combo (Max of Cor Healthy Cor Cancer then Squared for each Gene)
                                maxCor <- pmax(mcols(linksx)$CorrelationHealthy,mcols(linksx)$CorrelationDisease,na.rm=TRUE)
                                linkageScore <- split(maxCor^2, mcols(linksx)$gene_name) %>% lapply(., sum) %>% unlist(use.names=TRUE)
				#Return Summary Data Frame
                                data.frame(TF = tfs[x], Gene = names(nSubpop), N = nSubpop, P = nSubpop / ncol(sigMatP2Gx), linkageScore = linkageScore[names(nSubpop)])
			},error=function(e){cat(sprintf("INFO : Invalid TF : [ %s ]\n",tfs[x]))})
})

t2gDF=do.call(rbind,t2gDF)
saveRDS(t2gDF,"t2gDF.rds")
