CHAIN="/home/ye/Work/BioAligment/SNP/Shi/chainTO/GRCh37_to_GRCh38.chain.gz"
DIR="/home/ye/Work/BioAligment/SNP/Shi/svm_explain/Disease"
PATTERN="coord.bed"
for bed in `find $DIR -name $PATTERN`
do
	path=`dirname $bed`
	output=${path}/CRCH38.bed
	unmap=${path}/umapped.bed
	/home/ye/anaconda3/envs/scatac/bin/liftOver $bed $CHAIN $output $unmap
done
