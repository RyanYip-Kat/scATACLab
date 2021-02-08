cellranger="/home/ye/Software/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger"
rsrc="/home/ye/Work/R/scATAC/ArchR/20201010/src-developing/export_seurat_cells.R"
rscript="/home/ye/anaconda3/envs/r-latest/bin/Rscript"

seurat_list=$1
out=$2
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
	origin_h5=${arr[2]}
        output=${out}/$sample
        $rscript $rsrc --seurat $seurat --outdir $output --rename
	$cellranger reanalyze --id $sample --matrix $origin_h5 --barcodes $output/barcode.csv   
done

