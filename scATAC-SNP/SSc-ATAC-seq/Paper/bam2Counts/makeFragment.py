import numpy as np
import pandas as pd
import argparse
import os
import subprocess

#################### configure
samtools="/home/ye/anaconda3/envs/BulkBio/bin/samtools"
macs2="/home/ye/anaconda3/envs/BulkBio/bin/macs2"
sinto="/home/ye/anaconda3/envs/BulkBio/bin/sinto"


#################### functions
def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def sortBam(inBam,outBam,threads=8):
    """sort and index bam file."""
    print("INFO : Sort bam ...")
    cmd="{} sort -@ {} {} -o {}".format(samtools,str(threads),inBam,outBam)
    submitter(cmd)
    print("INFO : Build bam index ...")
    cmd="{} index -@ {} {}".format(samtools,str(threads),outBam)
    submitter(cmd)
    print("INFO : Done!")

def makeFragment(inBam,outFile="fragments.tsv",threads=8):
    print("INFO : Make Fragment file ...")
    cmd="{} fragments -b {} -p {} -f {}".format(sinto,inBam,str(threads),outFile) + " " + "--barcode_regex" + " " + '"[^:]*"'
    print("INFO : ---  [ {} ]".format(cmd))
    submitter(cmd)
    print("INFO : Done!")


def callpeak(inBam,genome_size=2.7e9,name="BC",outDir="callpeak",shift=None,extsize=None,method="q"):
    """
    inBam : input bam file
    genome_size :  genome_size (humen or mouse)
    name : prefix name for *_summits.bed
    outDir : output path
    """
    makedir(outDir)
    print("INFO : macs2 call peak ...")
    cmd="{} callpeak -g {} --name {} --treatment {} --outdir {} --format BAM --nomodel --call-summits --nolambda --keep-dup all".format(
            macs2,str(genome_size),str(name),inBam,outDir)

    if shift is not None and  extsize is not None:
        cmd=cmd + " " + " --shift {} --extsize {}".format(str(shift),str(extsize))

    if method=="q":
        cmd=cmd + " " + "-p" + " " + str(cutoff)
    else:
        cmd=cmd + " " + "-q" + " " + str(cutoff)
    print("INFO : ---  [ {} ]".format(cmd))
    submitter(cmd)
    print("INFO : Done!")


def get_args():
    parser = argparse.ArgumentParser("get counts from bam file")
    parser.add_argument('--bams',type=str,nargs="+",default=None, help='path of  bam files')
    parser.add_argument('--jobs',type=int,default=12, help='the number threads')
    parser.add_argument('--outdir',type=str,default="./Results", help='the path to save result')
    args = parser.parse_args()
    return args

if __name__=="__main__":
   
    args=get_args()
    bams=args.bams
    outDir=args.outdir
    threads=args.jobs
    makedir(outDir)
    
    print("INFO : Start work ...")
    for bam in bams:
        try:
            name=os.path.basename(bam).split(".")[0]
            outBam=os.path.join(os.path.dirname(bam),name+".sorted.bam")
            if not os.path.exists(outBam):
                sortBam(bam,outBam,threads) # sort bam file
            outFrag=os.path.join(outDir,name+"_fragments.tsv")
            makeFragment(outBam,outFrag,threads)
        except:
            print("Please check your function or bam file!!!")

    print("INFO : All Done!")









