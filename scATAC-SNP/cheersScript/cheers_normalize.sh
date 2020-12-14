src="/home/ye/Work/BioAligment/SNP/CHEERS/CHEERS_normalize.py"
python="/home/ye/anaconda3/envs/ldsc/bin/python"

#FeatureCount_DIR="/home/ye/Work/BioAligment/SNP/Shi/featureCounts/label_fine"
#OUT="./LD_Cheers_Normalize/"
FeatureCount_DIR=$1
OUT=$2
if [ ! -d $OUT ];
then
        mkdir -p $OUT
fi

$python $src --input $FeatureCount_DIR  --outdir $OUT --prefix "cheers" 
