#SNP_PATH="/home/ye/Work/BioAligment/SNP/snp-list"
SNP_PATH="./SNP"
CHAIN="/home/ye/Work/BioAligment/SNP/Shi/chainTO/hg19ToHg38.over.chain"
for snp in `ls $SNP_PATH`
do
	SNP=${SNP_PATH}/${snp}
	echo $SNP
	/home/ye/anaconda3/envs/scatac/bin/Rscript /home/ye/Work/BioAligment/SNP/Shi/CHEERS/snp_ld_cheersBlock_coord_transform.R --chain $CHAIN --snp $SNP --outdir ./snpHg38
done
