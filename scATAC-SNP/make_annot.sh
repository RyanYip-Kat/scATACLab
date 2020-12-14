src="/home/ye/Work/BioAligment/SNP/ldsc/make_annot.py"
python="/home/ye/anaconda3/envs/ldsc/bin/python"

OUT="LD_Annot"
if [ ! -d $OUT ];
then
	mkdir $OUT
fi


NAME="GDT_Cell"
BED_DIR="/home/ye/Work/BioAligment/SNP/featureCount"
BIM_DIR="/home/ye/Data/SNP_List/LDSCfile/1000G_Phase3_EAS_plinkfiles"

BED_FILES=`ls $BED_DIR | grep bed`
for chr in {1..22}
do      
	echo "---------"
	echo "### Make Annot File with Chrom :  $chr"
	for bed in $BED_FILES
	do
		outfile=$OUT/$bed
		if [ ! -d $outfile ];
		then
			mkdir $outfile
		fi

		bim_file=${BIM_DIR}/1000G.EAS.QC.${chr}.bim
		bed_file=${BED_DIR}/$bed
		annot=$outfile/${NAME}.${chr}.annot.gz
		echo "***  With BED File : $bed ****"
		$python $src --bed-file $bed_file --bimfile $bim_file --annot-file $annot
	done
done


