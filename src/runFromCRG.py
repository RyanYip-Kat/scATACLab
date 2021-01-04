import subprocess
import os
import argparse

#######################
samtools="/home/ye/anaconda3/envs/SnapATAC/bin/samtools"
snaptools="/home/ye/anaconda3/envs/SnapATAC/bin/snaptools"


####################### functions 
def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

def mkdir(DirX):
    if not os.path.exists(DirX):
        os.makedirs(DirX)


def prepareBam(bam,outDir,rmTemp=False):
    mkdir(outDir)
    result={}

    print("INFO : extract the header file ...")
    header_sam=outDir+"/"+"possorted.header.sam"
    result["header_sam"]=header_sam
    if not os.path.exists(header_sam):
        cmd="{} view {} -H > {}".format(samtools,bam,header_sam)
        submitter(cmd)

    print("INFO : create a bam file with the barcode embedded into the read name ...")
    cmd="bash ./src/createBam.sh"
    snap_bam=outDir+"/"+"possorted.snap.bam"
    result["snap_bam"]=snap_bam
    if not os.path.exists(snap_bam):
        cmd+=" {} {} {} {}".format(samtools,header_sam,bam,snap_bam)
        submitter(cmd)

    print("INFO : sort bam ...")
    sort_bam=outDir+"/"+"possorted.snap.nsrt.bam"
    result["sort_bam"]=sort_bam
    if not os.path.exists(sort_bam):
        cmd="{} sort -n -@ 10 {} -o {}".format(samtools,snap_bam,sort_bam)
        submitter(cmd)
    
    # whether or not remove temp
    if rmTemp:  
        cmd="rm -r {}".format(header_sam)
        submitter(cmd)
        cmd="rm -r {}".format(snap_bam)
        submitter(cmd)

    return sort_bam

def prepareFragment(fragment,outDir):
    mkdir(outDir)
    result={}
    if fragment.endswith(".gz"):
        print("INFO : unzip file ...")
        cmd="gunzip {}".format(fragment)
        submitter(cmd)
        name=os.path.basename(fragment)
        fragment=os.path.join(os.path.dirname(fragment),name.replace(".gz",""))
    
    fragment_bed=outDir+"/"+"fragments.bed"
    fragment_bed_gzip=fragment_bed+".gz"
    
    if not os.path.exists(fragment_bed_gzip):
        if not os.path.exists(fragment_bed):
            cmd="sort -k4,4 {} > {}".format(fragment,fragment_bed)
            submitter(cmd)
            cmd="gzip {}".format(fragment_bed)
            submitter(cmd)
        else:
            cmd="gzip {}".format(fragment_bed)
            submitter(cmd)
    
    return fragment_bed_gzip

def snap_pre(inFile,
        genome="mm10",
        genome_size="mm10.chrom.sizes",
        outSnap="crg.snap"):

    cmd="{} snap-pre --input-file {} --output-snap {} --genome-name {} --genome-size {} \
            --min-mapq 30 --min-flen 50 --max-flen 1000 --keep-chrm False --keep-single False \
            --keep-secondary False  --overwrite True --max-num 20000 --min-cov 500  --verbose True".format(snaptools,
                    inFile,outSnap,genome,genome_size)
    if not os.path.exists(outSnap):
        print("INFO : get snap file ...")
        submitter(cmd)


def snap_addBmat(snap_file):
    cmd="{} snap-add-bmat --snap-file {} --bin-size-list 1000 5000 10000 --verbose True".format(snaptools,snap_file)
    submitter(cmd)
   

def get_args():
    parser = argparse.ArgumentParser(description='Program to handle cellranger atac output width sanptools')
    parser.add_argument('--input',type=str,default=None,help="cellranger atac output's bam file")
    parser.add_argument("--outdir",type=str,default="./Results",help="path to save results")

    parser.add_argument("--genome",type=str,default="mm10",choices=["mm10","hg19","hg38"],help="genome name")
    parser.add_argument("--genome_size",type=str,default="mm10.chrom.sizes",help="genome chrom size file")

    parser.add_argument("--isbam",action="store_true",default=False,help="whether is bam file(or fragment")
    parser.add_argument("--cleanTmp",action="store_true",default=False,help="whether reomve temp bam file")
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args=get_args()

    inFile=args.input
    outDir=args.outdir
    genome=args.genome
    genome_size=args.genome_size
    isBam=args.isbam
    cleanTmp=args.cleanTmp
    mkdir(outDir)
    
    print("INFO : Prepare ...")
    if isBam:
        preSnap=prepareBam(inFile,outDir,rmTemp=cleanTmp)
    else:
        preSnap=prepareFragment(inFile,outDir)

    print("INFO : snap-pre ...")
    outSnap=outDir+"/"+"crg.snap"
    snap_pre(preSnap,genome=genome,genome_size=genome_size,outSnap=outSnap)

    print("INFO : add Bmat ...")
    snap_addBmat(outSnap)
    print("INFO : Done!")
