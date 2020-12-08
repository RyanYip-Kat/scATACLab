library(Seurat)
library(ggplot2)
library(argparse)
library(ggpubr)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="output",
                    help="the path to save result")

parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="the path to save result")

#parser$add_argument("--genes",
#		    nargs="+",
#                    type="character",
#                    default=NULL)

parser$add_argument("--groupby",
		    type="character",
                    default="condition")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

seurat=readRDS(args$seurat)
#genes=args$genes
genes=rownames(seurat)
metadata=seurat@meta.data
if(!args$groupby%in%colnames(metadata)){
	stop("Invalid group use !!!")
}

print(paste0("Violin Plot with : ",length(genes)," genes"))

DATA=GetAssayData(seurat,"data")
Idents(seurat)=metadata[[args$groupby]]
groups=as.character(unique(Idents(seurat)))
cb=combn(groups,2)

Gene_Pdata=list()
k=1
for(gene in genes){
	data=list()
	for(g in groups){
		cell=WhichCells(seurat,idents=g)
	        df=DATA[gene,cell]
		data[[g]]=df
	}
	
	Names=c()
	ps=c()
	for(i in 1:ncol(cb)){
		x_name=cb[1,i]
		y_name=cb[2,i]
		x=data[[x_name]]
		y=data[[y_name]]
		name=paste0(x_name,"_",y_name)
		pvalue=wilcox.test(x,y)$p.value
		ps<-c(ps,pvalue)
		Names<-c(Names,name)
	}
	gene_p<-data.frame(t(ps))
	colnames(gene_p)=Names
	rownames(gene_p)=gene
	Gene_Pdata[[k]]=gene_p
	k=k+1
}

df=do.call(rbind,Gene_Pdata)
write.table(df,file.path(args$outdir,"genes_pvalue_table.csv"),sep=",",quote=F)
