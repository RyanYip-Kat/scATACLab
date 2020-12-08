library(argparse)
library(stringr)
#library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the project path  of ArchR")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--column",
                    type="character",
                    default="the column in cellcoldata in archr as group")

parser$add_argument("--label_dict",
                    type="character",
                    default="the csv file of label dict")

parser$add_argument("--name",
                    type="character",
                    default="the name add  on cellcoldata in archr as group")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
projHeme=loadArchRProject(args$project)


metadata=as.data.frame(getCellColData(projHeme))
stopifnot(args$column%in%colnames(metadata))

dict=read.csv(args$label_dict,stringsAsFactors=FALSE)
stopifnot(ncol(dict)==2)
colnames(dict)=c("old_labels","new_labels")

target_label=metadata[[args$column]]
labelNew=dict$new_labels
labelOld=dict$old_labels

clusters=mapLabels(target_label, newLabels = labelNew, oldLabels = labelOld)
metadata$Labels=clusters
print(table(clusters))
print(table(target_label))

print(paste0("# Add name :",args$name," on cellcolmetadata"))
projHeme=addCellColData(ArchRProj=projHeme,data=clusters,name=args$name,cells=rownames(metadata))
saveArchRProject(ArchRProj = projHeme, outputDirectory =args$outdir , load = TRUE)
