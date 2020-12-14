rscript="/home/ye/anaconda3/envs/scatac/bin/Rscript"
snp_gr="/home/ye/Work/BioAligment/SNP/Shi/master_SNP/hg38_Diy-GWAS_SNPS.gr.rds"
conn_code="/home/ye/Work/BioAligment/SNP/Shi/scripts/PeakSet_Run_Cicero.R"   # --project --outdir --binarize
gwas_snp_co="/home/ye/Work/BioAligment/SNP/Shi/scripts/ChromVAR_For_GWAS_w_CoAccessbility_v2.R" # --se --conn --SNP --cutoff --outdir 
plot_code="/home/ye/Work/BioAligment/SNP/Shi/scripts/deviationScoresHeatmap.R" # --dev --outdir --groupby 

project=$1
outdir=$2

cutoff=0.35
groupby="label_fine"
echo "### Run Cicero"
$rscript $conn_code --project $project --outdir $outdir --binarize

echo "### Run GWAS _CoAccessbility"
$rscript $gwas_snp_co --se ${outdir}/scATAC-Summarized-Experiment.rds --conn ${outdir}/Peaks-Co-Accessibility.rds --SNP $snp_gr --cutoff $cutoff --outdir $outdir

echo "### Plot Deviatation Score's heatmap"
$rscript $plot_code --dev ${outdir}/GWAS-Co-Accessibility-chromVAR-Summarized-Experiment.rds --outdir $outdir --groupby $groupby


