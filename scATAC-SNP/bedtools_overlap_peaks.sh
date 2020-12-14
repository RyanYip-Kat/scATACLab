usage()
{
        printf "Usage: overlap peak [options]\n\n"
        printf "Options:\n"
        printf "\t-o OUTDIR , --outdir OUTDIR\n"
        printf "\t\tThe path for saving result\n"
        printf "\t-a QUERY\n"
	printf "\t\tThe SNP query bed file for bedtools overlap\n"
        printf "\t-b REFERENCE\n"
	printf "\t\tThe refernce bed(Cluster) file using for bedtools overlap\n"
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
                -b )
                        REF=$VALUE
                        shift 2
                        ;;
                -a )
                        QUERY=$VALUE
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

out=$OUTDIR/overlap.bed
time $BEDTOOLS intersect -a $QUERY -b $REF > $out


