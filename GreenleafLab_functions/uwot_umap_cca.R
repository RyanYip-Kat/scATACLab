library(argparse)
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--outdir",
                    type="character",
                    default="./output")


parser$add_argument("--cca",
                    type="character",
                    default=NULL,
                    help="cca aligned from scATAC genescore and scRNA aligned")


args <- parser$parse_args()
if(!dir.exists(args$outdir)){
        dir.create(args$outdir,recursive=TRUE)
}

set.seed(1)
print("### Loading CCA Aligned ")
alignedCCA=readRDS(args$cca)

print("### Run UMAP")
umap <- uwot::umap(
    alignedCCA, 
    n_neighbors = 50, 
    min_dist = 0.5, 
    metric = "euclidean", 
    n_threads = 5, 
    verbose = TRUE, 
    ret_model = FALSE)

#Plot DF
saveRDS(umap,file.path(args$outdir,"alignedCCA_umap.rds"))
