library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(argparse)
library(stringr)
library(future)
library(cicero)
#############################
parser <- ArgumentParser(description='Program For Finding co-accessible networks with Cicero')
parser$add_argument("--seurat",
                    type="character",
                    default=NULL,
                    help="Seurat(Signac) Object rds file")


parser$add_argument("--outdir",
                    type="character",
                    default="./Results")


args <- parser$parse_args()
############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}


############################## Configure
plan("multiprocess", workers = 16)
outDir=args$outdir
makedir(outDir)


############################## Loading dataset and convert into Cicero Object
message("INFO : Loading dataset ...")
seurat=readRDS(args$seurat)

message("INFO : as.cell_data_set  ...")
seurat.cds <- as.cell_data_set(x = seurat,assay="ATAC")

message("INFO : Convert into Cicero Object ...")
seurat.cicero <- make_cicero_cds(seurat.cds, reduced_coordinates = reducedDims(seurat.cds)$UMAP)

message("INFO : get the chromosome sizes from the Seurat object ...")
genome <- seqlengths(seurat)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
#genome <- genome[1]

message("INFO : convert chromosome sizes to a dataframe ...")
genome.df <- data.frame("chr" = names(genome), "length" = genome)

message("INFO : run cicero ...")
conns <- run_cicero(seurat.cicero, genomic_coords = genome.df, sample_num = 100)
saveRDS(conns,file.path(outDir,"connections.rds"))

message("INFO : Find cis-co-accessible networks (CCANs) ...")
ccans <- generate_ccans(conns)
saveRDS(ccans,file.path(outDir,"ccans.rds"))

message("INFO : Add links to a Seurat object ...")
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(seurat) <- links

message("INFO : Save Object ...")
saveRDS(seurat,file.path(outDir,"seurat-cis-links.rds"))
