#!/bin/bash

SCRIPT_PATH="/home/ye/Work/BioAligment/PacificBiosciences/scisoseq-counts/bin"
CELLRANGER_PY="/home/ye/Software/cellranger-3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin/python"
CELLRANGER_PATH="/home/ye/Software/cellranger-3.1.0"

usage()
{
        printf "Usage: Isoseq3 pipeline [options]\n\n"
        printf "Options:\n"
        printf "\t--output-dir OUTPUT-DIR\n"
        printf "\t\tthe output dir to save result.\n"
        printf "\t--cellranger-matrix\n"
        printf "\t\tthe matrix.h5  from cellranger\n"
	printf "\t--target-matrix\n"
	printf "\t\tthe *.h5 need to be converted\n"
        printf "\t-h, --help\n"
        printf "\t\tShow this help message and exit.\n"
}

while [ "$1" != "" ]; do
        PARAM=`echo $1 | awk -F' ' '{print $1}'`
        VALUE=`echo $2 | awk -F' ' '{print $1}'`
        case $PARAM in
                -h | --help)
                        usage
                        exit
                        ;;
                --cellranger-matrix)
                        CRH5=$VALUE
                        shift 2
                        ;;
		--target-matrix)
			FROM_H5=$VALUE
			shift 2
			;;
                --output-dir)
                        OUT=$VALUE
                        shift 2
                        ;;
                *)
                        echo "Error: unknown parameter \"$PARAM\""
                        exit 1
                        ;;
        esac
done


if [ ! -d $OUT ];
then
	mkdir -p $OUT
fi

echo "Convert matrix to CellRanger h5 format."
$CELLRANGER_PY $SCRIPT_PATH/scMatrix2CellRangerH5.py -m $FROM_H5  -c $CRH5 -o $OUT/filtered_feature_bc_matrix_gold.h5 --cellranger_path $CELLRANGER_PATH


