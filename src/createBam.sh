samtools=$1
header_sam=$2
possorted_bam=$3
outBam=$4
cat <( cat $header_sam ) <( $samtools view $possorted_bam  | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) | $samtools view -bS - > $outBam
