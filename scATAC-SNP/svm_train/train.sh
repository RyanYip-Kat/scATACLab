path=$1  # must absolute path
#out=$2

for name in `ls $path`
do
	pos=${path}/${name}/ctcfpos.fa
	neg=${path}/${name}/ctcfneg.fa
	prefix=$name
	outdir=${path}/${name}
	if [ ! -d $outdir ]
	then
           mkdir -p $outdir
        fi
        
	cd $outdir
	echo "### Train lSK SVM with :  $name"
	/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmtrain $pos $neg $prefix -x 5 -i 3 -T 16 -r 9527 -m 5000 
	#/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmpredict $pos ${prefix}.model.txt ${name}.pos.txt -T 16 
	#/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmpredict $neg ${prefix}.model.txt ${name}.neg.txt -T 16
	cd -
done
	
