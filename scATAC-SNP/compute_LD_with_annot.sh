src="/home/ye/Work/BioAligment/SNP/ldsc/ldsc.py"
python="/home/ye/anaconda3/envs/ldsc/bin/python"

#OUT="annot_files"
#if [ ! -d $OUT ];
#then
#	mkdir $OUT
#fi


NAME="GDT_Cell"
Annot_DIR="/home/ye/Work/BioAligment/SNP/LD_Annot"
BFILE_DIR="/home/ye/Data/SNP_List/LDSCfile/1000G_Phase3_EAS_plinkfiles"
HM_SNP="/home/ye/Data/SNP_List/LDSCfile/hapmap3_snps"

for chr in {1..22}
do      
	echo "------------- Chrom : chr_${chr}"
	for dir in `ls $Annot_DIR`
	do
		annot=$Annot_DIR/$dir/${NAME}.${chr}.annot.gz
		snp=$HM_SNP/hm.${chr}.snp
		bfile=$BFILE_DIR/1000G.EAS.QC.${chr}
		outfile=$Annot_DIR/$dir/${NAME}.${chr}

		echo "Compute LD Score in Cell : $dir"
		echo "      With SNP : $snp"
		echo "      With bfile : $bfile"
		echo "      With Annot : $annot"
		echo "      Output : $outfile"
		$python $src --l2 --bfile $bfile --ld-wind-cm 1 --annot $annot  --thin-annot --out $outfile --print-snps $snp
	done
done


