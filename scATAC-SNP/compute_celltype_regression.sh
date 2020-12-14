src="/home/ye/Work/BioAligment/SNP/ldsc/ldsc.py"
python="/home/ye/anaconda3/envs/ldsc/bin/python"


NAME="GDT_Cell"
baseline="/home/ye/Data/SNP_List/LDSCfile/1000G_Phase3_EAS_baselineLD_v2"
#sumstats="/home/ye/Work/BioAligment/SNP/snp-sumstat/PASS_Type_2_Diabetes.sumstats.gz"
sumstats="/home/ye/Work/BioAligment/SNP/cooper_2012.sumstats.gz"
ldct="/home/ye/Work/BioAligment/SNP/GDT_Cell.ldcts"
weight="/home/ye/Work/BioAligment/SNP/ldsc-seg/weights_hm3_no_hla"

$python $src --h2-cts $sumstats \
	--ref-ld-chr $baseline/baselineLD. \
	--out $NAME \
	--ref-ld-chr-cts $ldct  \
	--w-ld-chr $weight/weights. 
	#--n-blocks 1000 
