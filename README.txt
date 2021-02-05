####### first
#   download http://archive.gersteinlab.org/proj/PeakSeq/Mappability_Map/Code 
#   download mm10_cisbp_fimo.zip(or huamn) reference database
#   install bagfoot package (  install.package("bagfoot_0.9.7.07.tar.gz"))

####### prepare data
# merge bam file : bash scripts/mergeBamFiles.sh 
# generate *_withMap.txt file and nuccode_*_6mer_chr*.dat  files  : Rscripts scripts/bagfootBam.R 
# generate hot2pot csv file : Rscripts scripts/ArchR2hotspot.R 
# generate bagfoot result : Rscripts scripts/bagfootRun.R



########################
reference link : https://sourceforge.net/projects/bagfootr/files/
