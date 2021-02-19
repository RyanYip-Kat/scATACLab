library(argparse)
library(stringr)
library(ggrepel)
library(ArchR)
library(Matrix)
library(ggplot2)

source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Trajectory plot')
parser$add_argument("--project",
                    type="character",
                    default="the path of project saved")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")


parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

parser$add_argument("--groupby",
		    type="character",
		    default=NULL,
		    help="column in metadata as group")

args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}


options(ggrepel.max.overlaps = Inf)
print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 


print("### Loading ArrowFiles")
projHeme5<-loadArchRProject(args$project)

message("INFO : get Motif GroupSE")
seGroupMotif <- getGroupSE(ArchRProj = projHeme5, useMatrix = "MotifMatrix", groupBy = args$groupby)
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

message("INFO : Identify Correlated TF Motifs and TF Gene Score/Expression")
corGSM_MM <- correlateMatrices(
    ArchRProj = projHeme5,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"

Text=with(corGSM_MM,ifelse(TFRegulator=="YES",paste(GeneScoreMatrix_seqnames,MotifMatrix_seqnames,GeneScoreMatrix_matchName,sep=":"),""))
corGSM_MM$Text=Text
corGSM_MM <- as.data.frame(corGSM_MM)
#p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
#  geom_point() +
#  geom_text_repel(data=corGSM_MM,aes(x = cor,y=maxDelta,label=Text),size = 5,box.padding = unit(0.35, "lines"),
#                  point.padding = unit(0.5, "lines"),
#                  segment.color = "grey50",
#                  show.legend = TRUE,
#                  colour = "black")+
#  theme_ArchR() +
#  geom_vline(xintercept = 0, lty = "dashed") +
#  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
#  xlab("Correlation To Gene Score") +
#  ylab("Max TF Motif Delta") +
#  scale_y_continuous(
#    expand = c(0,0),
#    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
#  )


p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator,label=Text)) +
  geom_point() +
  geom_text_repel(size = 5,box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = "grey50",
                   show.legend = TRUE,
                   colour = "black")+
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )


ggsave(file.path(args$outdir,"PosTFRegulator.pdf"),plot=p,width=16,height=12)
write.table(corGSM_MM,file.path(args$outdir,"corGSM_MM.csv"),sep=",",quote=F,row.names=F)
