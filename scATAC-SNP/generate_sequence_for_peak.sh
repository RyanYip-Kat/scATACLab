usage()
{
        printf "Usage: biomex [options]\n\n"
        printf "Options:\n"
        printf "\t-o OUTDIR , --outdir OUTDIR\n"
        printf "\t\tThe path for saving result\n"
        printf "\t-fa FASTA, --fasta FASTA\n"
	printf "\t\tThe genome fasta file for generate sequence\n"
        printf "\t-b BED,--bed BED\n"
        printf "\t\tThe bed file using for generating\n"
        printf "\t-h, --help\n"
        printf "\t\tShow this help message and exit.\n"
}

# Get the parameters selected by the user
while [ "$1" != "" ]; do
        PARAM=`echo $1 | awk -F' ' '{print $1}'`
        VALUE=`echo $2 | awk -F' ' '{print $1}'`
        case $PARAM in
                -h | --help)
                        usage
                        exit
                        ;;
                -b | --bed)
                        BED=$VALUE
                        shift 2
                        ;;
                -fa | --fasta)
                        FASTA=$VALUE
                        shift 2
                        ;;
                -o | --outdir)
                        OUTDIR=$VALUE
                        shift 2
                        ;;
                *)
                        echo "Error: unknown parameter \"$PARAM\""
                        exit 1
                        ;;
        esac
done

BEDTOOLS="/home/ye/anaconda3/envs/BulkBio/bin/bedtools"
if [ ! -d $OUTDIR ];
then
	mkdir -p $OUTDIR
fi

outFA=$OUTDIR/sequence.fa
echo "### The BED File From : $BED"
echo "### The FASTA GENOME From : $FASTA"
echo "### The Result Save into :  $outFA"

time $BEDTOOLS getfasta -fi $FASTA -bed $BED -fo $outFA


