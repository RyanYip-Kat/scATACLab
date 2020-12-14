nfolds=5
workers=24
out="./gkmModels"
for bed in `cat ./reproduciblePeaks_15-28.txt`
do
	file=`basename $bed`
	name=${file%.*}
	printf "INFO : $name : $bed\n"
	python main.py --bed $bed --nfold $nfolds --outdir ${out}/${name} --workers $workers 
done
