project=$1
outdir=$2
/home/ye/anaconda3/envs/scatac/bin/Rscript /home/ye/Work/R/scATAC/ArchR/multiome-atac-gex/getDoubletScore.R --project $project --outdir $outdir
