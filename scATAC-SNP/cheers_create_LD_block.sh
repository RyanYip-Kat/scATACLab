src="/home/ye/Work/BioAligment/SNP/CHEERS/create_LD_blocks.py"
python="/home/ye/anaconda3/envs/ldsc/bin/python"

LD_DIR="/home/ye/Data/SNP_List/cheers"
SNP_DIR="/home/ye/Work/BioAligment/SNP/imd/snp_list"

OUT="LD_block"
if [ ! -d $OUT ];
then
        mkdir $OUT
fi

for list in `ls $SNP_DIR`
do
	name=${list%.*}
	snp_list=$SNP_DIR/$list
	outdir=$OUT/$name
	if [ ! -d $outdir ];
	then
		mkdir $outdir
	fi

        echo "##### Create LD Block :"
	echo "      SNP LIST : $snp_list"
	echo "      OUTPUT : $outdir"
	echo "      LD DIR : $LD_DIR"
	time $python $src $snp_list $outdir $LD_DIR 
done	
