write_src="/home/ye/Work/R/scATAC/ArchR/utils/write10xCounts_h5_from_seurat.R"
transform_src="/home/ye/Work/R/scATAC/ArchR/utils/transform_h5.sh"
rscript="/home/ye/anaconda3/envs/r-latest/bin/Rscript"

seurat_list=$1
out=$2
CRH5=$3
if [ ! -d $out ]
then
	mkdir -p $out
fi

cat $seurat_list | while read id
do
	arr=(${id})
	sample=${arr[0]}
	echo " Convert Sample : $sample"
	seurat=${arr[1]}
	output=${out}/$sample
	$rscript $write_src --seurat $seurat --outdir $output --rename 
	bash $transform_src --target-matrix ${output}/matrix --cellranger-matrix $CRH5 --output-dir $output
done

