snap=$1
/home/ye/anaconda3/envs/SnapATAC/bin/snaptools snap-add-bmat  \
	--snap-file=$snap  \
	--bin-size-list 1000 5000 10000  \
	--verbose=True
