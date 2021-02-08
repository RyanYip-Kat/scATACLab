config=$1
output=$2

python="/home/ye/anaconda3/envs/BulkBio/bin/python"
cat $config  | while read id
do
	arr=(${id})
	sample=${arr[0]}
	path=${arr[1]}
	echo $sample  "######"  $path
	output_path=$output/$sample
	if [ ! -d $output_path ]
	then
		mkdir -p $output_path
        fi

	$python /home/ye/Work/Python/SingleCell/Project/Paper-scATAC-scRNA/src/atac/clean_barcode_multiplets_1.1.py --input_path $path --output_path $output_path
done
