library(argparse)
library(stringr)
library(ArchR)
library(Matrix)
library(ggplot2)
library(Cairo)
source("/home/ye/Work/R/scATAC/ArchR/plotter/plotDF.R")
#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
		    help="the project path  of ArchR")


parser$add_argument("--outdir",
                    type="character",
                    default="ArchR_result")

parser$add_argument("--group",
                    type="character",
                    default="Clusters")

parser$add_argument("--width",
                    type="integer",
                    default=12,
                    help="the width of plot")

parser$add_argument("--height",
                    type="integer",
                    default=12,
                    help="the height of plot")



args <- parser$parse_args()

if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}
projHeme=loadArchRProject(args$project)


width=args$width
height=args$height

print(paste0("# The plot size is : [ ",width,height," ]"))
markerTest <- getMarkerFeatures(
  ArchRProj = projHeme,
  useMatrix = "PeakMatrix",
  groupBy = args$group,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)

print("# motif up ")
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

print("# motif down")
motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

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

MyplotPDF(ggUp, ggDo, name = paste0(args$group,"-Markers-Motifs-Enriched"), width = width, height = height,outpath=args$outdir, addDOC = FALSE)


enrichMotifs <- peakAnnoEnrichment(
    seMarker =markerTest,
    ArchRProj = projHeme,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )


heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
MyplotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = width, height = height, outpath=args$outdir, addDOC = FALSE)
