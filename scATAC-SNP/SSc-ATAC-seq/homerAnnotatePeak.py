import subprocess
import os
import sys
import argparse

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='A program to AnnotatePeak bed file with annotatePeaks.pl in homer')
    parser.add_argument('bed')
    parser.add_argument('prefix')
    parser.add_argument('outdir')
    args = parser.parse_args()

    outdir=args.outdir
    Bed=args.bed
    prefix=args.prefix

    Annotate="/home/ye/anaconda3/envs/BulkBio/bin/annotatePeaks.pl"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print("INFO : Homer annotation Peaks")
    AnnoF=outdir + "/" + prefix + ".bed.ano.txt"
    Comd='{} {} hg38 > {}'.format(Annotate,Bed,AnnoF)
    submitter(Comd)
