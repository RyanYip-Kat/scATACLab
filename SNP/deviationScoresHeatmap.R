library(chromVAR)
library(dplyr)
library(SummarizedExperiment)
library(magrittr)
library(argparse)
library(stringr)
library(pheatmap)
parser <- ArgumentParser(description='Heatmap for deviationScores from GWAS-Co-Accessibility-chromVAR-Summarized-Experiment')
parser$add_argument("--dev",
                    type="character",
                    default=NULL,
                    help="chromVARDeviations from computeDeviation")

parser$add_argument("--outdir",
                    type="character",
                    default="./results")

parser$add_argument("--groupby",
                    type="character",
                    default="label_fine",
		    help="The Column in colData for heatmap plotting")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

message("INFO : Loading chromVARDeviations")
dev=readRDS(args$dev)
dev_score=deviationScores(dev)
metadata=as.data.frame(colData(dev))
Clusters=metadata[[args$groupby]]

score_name=rownames(dev_score)
score_list=list()
for(name in score_name){
	cat(sprintf("Get :%s median deviation score\n",name))
	score=dev_score[name,]
	X=data.frame("Score"=score,"Cluster"=Clusters)
	z=as.data.frame(X %>% group_by(Cluster) %>% summarise(v = median(Score,na.rm=TRUE)))
	rownames(z)=z$Cluster
	z=subset(z,select="v")
	colnames(z)=name
	score_list[[name]]=z
}

scores=do.call(cbind,score_list)

message("INFO : Save median deviatation score")
saveRDS(scores,file.path(args$outdir,paste0("Median-DeviationScores-Across-",args$groupby,".rds")))
scores[is.na(scores)]=0.0

pdf(file.path(args$outdir,paste0("Median-DeviationScores-Across-",args$groupby,".pdf")),width=12,height=16)
pheatmap(t(scores),
	 scale="column",
	 cluster_cols=TRUE,
	 cluster_rows = TRUE, 
	 fontsize = 10,
	 fontsize_row = 10,
	 fontsize_col=10
	 )

dev.off()

	



