src="/home/ye/Work/BioAligment/SNP/CHEERS/CHEERS_computeEnrichment.py"
python="/home/ye/anaconda3/envs/ldsc/bin/python"

CHEERS_normalize="/home/ye/Work/BioAligment/SNP/imd/LD_Cheers_Normalize/IMD_counts_normToMax_quantileNorm_euclideanNorm.txt" # from cheers_normalize
SNP_DIR="/home/ye/Work/BioAligment/SNP/imd/LD_block"  # from create_LD_block.py

OUT="./LD_Cheers_Enrichment/"
if [ ! -d $OUT ];
then
        mkdir $OUT
fi

for name in `ls $SNP_DIR`
do
        LD=$SNP_DIR/$name
	echo "### Compute CHEERS Enrichment"
	echo "    TRAIT  NAME : $name"
	echo "    LD TRAIT : $LD"
	$python $src --input $CHEERS_normalize  --ld $LD --outdir $OUT --trait $name
done
