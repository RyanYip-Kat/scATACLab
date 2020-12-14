BED=$1
CHAIN=$2
PREFIX=$3
OUT=$4
if [ ! -d $OUT ]
then
	mkdir -p $OUT
	
fi

filename=${PREFIX}.bed
unmap=${PREFIX}.unmap.bed
/home/ye/anaconda3/envs/scatac/bin/liftOver $BED $CHAIN $filename $unmap


