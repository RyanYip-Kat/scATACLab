library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
library(ggplot2)
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="the project path  of ArchR")

parser$add_argument("--logfc",
                    type="double",
                    default=1.0,
                    help="logfc cutoff for getMarkers")

parser$add_argument("--seMarker",
                    type="character",
                    default=NULL,
                    help="SummarizedExperiment object returned by markerFeatures")

parser$add_argument("--groupby",
                    type="character",
                    default="Clusters",
                    help="the column in cellcoldata in archr as group")

parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--column",
                    type="character",
                    default="Clusters",
                    help="the column in cellcol to be exported")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset of column to be extracted")


args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
message("INFO : Loading dataset")
logfc=abs(args$logfc)
projHeme=loadArchRProject(args$project)

metadata=as.data.frame(getCellColData(projHeme))
stopifnot(args$groupby%in%colnames(metadata))
stopifnot(args$useMatrix%in%getAvailableMatrices(projHeme))

cells=getCellNames(projHeme)
if(!is.null(args$column) & !is.null(args$subset)){
        stopifnot(args$column%in%colnames(metadata))
        target=metadata[[args$column]]

        stopifnot(args$subset%in%unique(target))
        idxPass= which(target%in%args$subset)
        cellPass=cells[idxPass]
        projHeme=subsetCells(projHeme,cellNames=cellPass)
}

if(is.null(args$seMarker)){
	seMarker <- getMarkerFeatures(
				ArchRProj = projHeme,
                                useMatrix = "PeakMatrix",
                                groupBy = args$groupby,
                                bias = c("TSSEnrichment", "log10(nFrags)"),
                                testMethod = "wilcoxon")

}else{
	seMarker=readRDS(args$seMarker)
}
message("INFO : Up  Enrichment ")
motifsUp <- peakAnnoEnrichment(
    seMarker = seMarker,
    ArchRProj = projHeme,
    peakAnnotation = "Motif",
    cutOff = paste0("Pval <= 0.05 & Log2FC >= ",logfc),
  )
saveRDS(motifsUp,file.path(args$outdir,"motifsUp.rds"))
message("INFO : Down  Enrichment ")
motifsDo <- peakAnnoEnrichment(
    seMarker = seMarker,
    ArchRProj = projHeme,
    peakAnnotation = "Motif",
    cutOff = paste0("Pval <= 0.05 & Log2FC <= ",-logfc),
  )

saveRDS(motifsDo,file.path(args$outdir,"motifsDo.rds"))
######################
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

write.table(df,file.path(args$outdir,"motifsUp.csv"),sep=",",row.names=F,quote=F)
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
write.table(df,file.path(args$outdir,"motifsDown.csv"),sep=",",row.names=F,quote=F)

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

MyplotPDF(ggUp,ggDo,name ="Markers-Motifs-Enriched", width = 16, height =12, outpath =args$outdir, addDOC = FALSE)

message("INFO : Done!")
