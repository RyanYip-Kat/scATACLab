library(ArchR)
library(ggplot2)
library(ggrepel)
library(stringr)
library(argparse)


#########################
softmax=function(V){
	expV=exp(V)
	return(expV/sum(expV))
}

makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

MinMaxScale<-function(x){
        stopifnot(is.vector(x))
        v=(x-min(x))/(max(x)-min(x)+1e-4)
        return(v)
}

rowZscores<-function (m = NULL, min = -2, max = 2, limit = FALSE)
{
    z <- sweep(m - Matrix::rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    if (limit) {
        z[z > max] <- max
        z[z < min] <- min
    }
    return(z)
}


DOCRPointPlot=function(project,corCutOff=0.45,resolution=1,nShowGenes=10,outDir="DORCPlot"){
	stopifnot(class(project)=="ArchRProject")
	makedir(outDir)
	genes=getFeatures(project,useMatrix="GeneIntegrationMatrix")
        p2g <- getPeak2GeneLinks(
				 ArchRProj = project,
				 corCutOff = corCutOff,
				 resolution = resolution,
				 returnLoops = FALSE)
        idxGenes=p2g$idxRNA
        idxGenesNum=table(idxGenes)
        idxGenesNum=idxGenesNum[order(idxGenesNum,decreasing=F)]
        df=data.frame("N"=as.integer(idxGenesNum),"R"=rank(idxGenesNum,ties.method="first"),"G"=names(idxGenesNum))
        df$genes=genes[df$G]

	nGenes=nShowGenes  #  number of gene to label
        label=df$genes
        label[1:(length(label)-nGenes)]=""
        df$Label=label

        noGrid <- theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), validate = TRUE)
	p=ggplot(data=df,aes(x=R,y=N))+
		geom_point()+
                geom_hline(yintercept = 10,linetype='dashed', color='black', size=2)+
                geom_vline(xintercept=min(df$R[df$N==10]),linetype='dashed', color='black', size=2)+
                xlim(c(0,max(df$R)+1000))+
                xlab("Ranked Sort Genes")+ylab("Number of Correlated Peaks")
	p=p+geom_text_repel(data=df,aes(x = R,y=N,label=Label),size = 5,box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = "grey50",
                  show.legend = TRUE,
                  colour = "black")+theme_bw()+noGrid

	ggsave(file.path(outDir,"DORCPointPlot.pdf"),plot=p,width=12,height=12)
}

getDORCXMatrix=function(project,colCellType="Clusters",dorcCutoff=10){
	#  project : ArchR Object with addPeak2Gene
	#  colCellType : column in metadata as celltype
	#  dorcCutoff :  the cutoff value for DORC

	stopifnot(class(project)=="ArchRProject")
	nFrags=project$nFrags
	peakset=getPeakSet(project)
	message("INFO :  get peak matrix ...")
	peakSE=getMatrixFromProject(project,useMatrix="PeakMatrix")
	peakCounts=assay(peakSE)  # peakCounts
	message("INFO : Normalize with unique fragments ...")
        peakCounts=peakCounts/nFrags*1e4

	geneMatrix=c("GeneIntegrationMatrix","GeneExpressionMatrix")[c("GeneIntegrationMatrix","GeneExpressionMatrix")%in%getAvailableMatrices(project)]
	genes=getFeatures(project,useMatrix=geneMatrix)
	message("INFO : get Peak2Gene Links ...")
        p2g <- getPeak2GeneLinks(
				 ArchRProj = project,
				 corCutOff = 0.45,
				 resolution = 1,
				 returnLoops = FALSE)
	
	idxGenes=p2g$idxRNA
        idxGenesNum=table(idxGenes)
        dorcName=as.integer(names(idxGenesNum[idxGenesNum > dorcCutoff]))
        dorcP2G=as.data.frame(p2g)[idxGenes%in%dorcName,]  #  number of gene correlated peaks

	message("INFO : get DORC DF ...")
	idxGenes=dorcP2G$idxRNA
        idxGenesNum=table(idxGenes)
        idxGenesNum=idxGenesNum[order(idxGenesNum,decreasing=F)]
        df=data.frame("N"=as.integer(idxGenesNum),"R"=rank(idxGenesNum,ties.method="first"),"G"=names(idxGenesNum))
        df$genes=genes[as.integer(df$G)]

	message("INFO : get DORC SCORE Matrix ...")
	dorcDF_list=list()
	for(i in 1:nrow(df)){
		G=df$G[i]
                gene=df$genes[i]
		if(i%%10==0){
			cat(sprintf("INFO : [ %d of %d ] --- [ DORC Gene idx [ %d of %s ]\n",i,nrow(df),G,gene))
		}
		idxG=dorcP2G$idxATAC[dorcP2G$idxRNA==G]
                peakG=peakCounts[idxG,]
                dorcG=colSums(peakG)
                dorcGDF=data.frame("DORC"=dorcG)
                colnames(dorcGDF)=gene
                dorcDF_list[[i]]=dorcGDF
	}
	dorcDF=do.call(cbind,dorcDF_list)  # cell x gene

	message("INFO : get Celltype DORC SCORE Matrix ...")
	celltypeDorcDF_list=list()
	metadata=getCellColData(project)
	CellType=metadata[[colCellType]]
	Cells=getCellNames(project)
	uCellType=unique(CellType)
	for(i in seq_along(uCellType)){
		celltype=uCellType[i]
		cat(sprintf("INFO : [ %s(%d) of %d ]\n",celltype,i,length(uCellType)))
		cell=Cells[CellType==celltype]
		cellMat=as.matrix(dorcDF[cell,])
		cellMatDF=data.frame("SCORE"=colMeans(cellMat)) # nGene x 1
		#cellMatDF=as.data.frame(apply(cellMatDF,2,softmax))
		colnames(cellMatDF)=celltype
		celltypeDorcDF_list[[i]]=cellMatDF
	}
	celltypeDorcDF=do.call(cbind,celltypeDorcDF_list)
	return(list("cellDORC"=dorcDF,"groupDORC"=celltypeDorcDF))
}
	
computeDORCKNN<-function (data,
			  k = 50,
			  knnIteration=500, 
			  includeSelf = FALSE,
			  overlapCutoff=0.8,
    ...)
{
    idx <- sample(seq_len(nrow(data)), knnIteration, replace = !nrow(data) >= knnIteration)
    query <- data[idx,]
    searchSelf <- FALSE

    require("nabor")
    if (searchSelf & !includeSelf) {
        knnIdx <- nabor::knn(data = data, query = query, k = k +
            1, ...)$nn.idx
        knnIdx <- knnIdx[, -1, drop = FALSE]
    }
    else {
        knnIdx <- nabor::knn(data = data, query = query, k = k,
            ...)$nn.idx
    }
    keepKnn <- ArchR:::determineOverlapCpp(knnIdx, floor(overlapCutoff * k))
    knnIdx <- knnIdx[keepKnn==0,]
    knnObj <- lapply(seq_len(nrow(knnIdx)), function(x){
    rownames(data)[knnIdx[x, ]]}) %>% SimpleList
    return(knnObj)
}



#projHeme=loadArchRProject("../multiome-atac-gex-test/Save-ProjHeme-P2G/")
#dorcList=getDORCXMatrix(projHeme)
#saveRDS(dorcList,"DORCScoreList.rds")

#x=dorcList[["CellTypeDORC"]]
#colData=data.frame("Clusters"=colnames(x))
#p=ArchR:::.ArchRHeatmap(mat=as.matrix(x),colData=colData,showRowDendrogram=TRUE,scale=TRUE,customRowLabel =c(1,3,5,7))
