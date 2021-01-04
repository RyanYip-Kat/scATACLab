snap=$1
peak=$2
/home/ye/anaconda3/envs/SnapATAC/bin/snaptools snap-add-pmat  \
	--snap-file=$snap  \
	--peak-file=$peak \
	--verbose=True
