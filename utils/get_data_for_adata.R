library(argparse)
library(stringr)
library(ArchR)
library(Matrix)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the ArchR project")

parser$add_argument("--outdir",
                    type="character",
                    default="output")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

print("# Loading ArchR project")
projHeme=loadArchRProject(args$project)

print("# Get gene score matrix")

features=getFeatures(projHeme,useMatrix="GeneScoreMatrix")
meta.data=as.data.frame(getCellColData(projHeme))
meta.data$cells=rownames(meta.data)
geneScore=getMatrixFromProject(projHeme,useMatrix="GeneScoreMatrix",verbose = TRUE)
rownames(geneScore)=features

print("# Convert into matrix")
score=assay(geneScore)

Genes<-data.frame(symbol=features)
score=as.data.frame(as.matrix(score))
counts<-cbind(Genes,score)

print("# Save")
write.table(counts,file.path(args$outdir,"cell_counts.txt"),sep="\t",row.names = FALSE,quote=F)
write.table(meta.data,file.path(args$outdir,"cell_meta.txt"),sep="\t",row.names = FALSE,quote=F)



