src="/home/ye/Work/BioAligment/SNP/CHEERS/CHEERS_normalize.py"
python="/home/ye/anaconda3/envs/ldsc/bin/python"

FeatureCount_DIR="/home/ye/Work/BioAligment/SNP/imd/featureCounts/label_main"

OUT="./LD_Cheers_Normalize/"
if [ ! -d $OUT ];
then
        mkdir $OUT
fi

$python $src --input $FeatureCount_DIR  --outdir $OUT --prefix "IMD" 
