import os
import sys
import pysam
import pandas as pd
import subprocess
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

ref_fasta="/home/ye/Data/10X/VDJ/ref/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa"
gen_dinucleotide_freqs="/home/ye/Work/BioAligment/SNP/alzheimers_parkinsons/svm_training/gc_dinuc_balanced/gen_dinucleotide_freqs.py"
Basedir="svm_train/positive/"
def get_seqs(bed, fasta):
    ref = pysam.FastaFile(ref_fasta)

    df_bed = pd.read_csv(bed, sep='\t', header=None)
    fa_file = open(fasta, 'w')

    counter = 0
    for index,row in df_bed.iterrows():
        seq=ref.fetch(row[0],int(row[1]),int(row[2]))
        fa_file.write('>' + str(counter) + '\n')
        fa_file.write(seq.upper() + '\n')
        counter += 1

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()


def train_preprocess(inputs):
    print('[TRAIN] Cluster: ' + str(inputs[0]) + '; Fold: ' + str(inputs[1]))
    basedir = os.path.join(Basedir,str(inputs[0]),'fold_' + str(inputs[1]), 'train')
    os.system('cut -f 1,2,3 ' + basedir + '/train.pos.bed | sed -e "s/$/\t1.0/" > ' + basedir + '/train.all.bed')
    #os.system("zcat " + basedir + "train.inputs.bed.gz | tail -n +2 | awk '$4 != 1.0' >> " + basedir + "train.all.bed")
    cmd="{} {} -b {} --ratio_neg_to_pos {} -o {} --ref_fasta {} --gc".format(
            "python",gen_dinucleotide_freqs,basedir + '/train.all.bed',1.0,
            basedir + '/train.gc.txt',ref_fasta)
    submitter(cmd)
    cmd="{} {} -b {} --ratio_neg_to_pos {} -o {} --ref_fasta {} --gc --dinuc_freq {}".format(
            "python",gen_dinucleotide_freqs,basedir + '/train.all.bed',1.0,
            basedir + '/train.final.bed',ref_fasta,basedir + '/train.gc.txt')
    
    submitter(cmd)
    pos_train = open(basedir + '/train.pos.bed')
    pos_train_len = len(pos_train.readlines())

    os.system('head -n ' + str(pos_train_len) + ' ' + basedir + '/train.final.bed > ' + basedir + '/train.final.pos.bed')
    os.system('tail -n +' + str(pos_train_len + 1) + ' ' + basedir + '/train.final.bed > ' + basedir + '/train.final.neg.bed')

    get_seqs(basedir + '/train.final.pos.bed', basedir + '/train.final.pos.fasta')
    get_seqs(basedir + '/train.final.neg.bed', basedir + '/train.final.neg.fasta')

    pos_train.close()

def test_preprocess(inputs):
    print('[TRAIN] Cluster: ' + str(inputs[0]) + '; Fold: ' + str(inputs[1]))
    basedir = os.path.join(Basedir,str(inputs[0]),'fold_' + str(inputs[1]), 'test')
    os.system('cut -f 1,2,3 ' + basedir + '/test.pos.bed | sed -e "s/$/\t1.0/" > ' + basedir + '/test.all.bed')
    #os.system("zcat " + basedir + "test.inputs.bed.gz | tail -n +2 | awk '$4 != 1.0' >> " + basedir + "test.all.bed")
    cmd="{} {} -b {} --ratio_neg_to_pos {} -o {} --ref_fasta {} --gc".format(
            "python",gen_dinucleotide_freqs,basedir + '/test.all.bed',1.0,
            basedir + '/test.gc.txt',ref_fasta)
    submitter(cmd)
    cmd="{} {} -b {} --ratio_neg_to_pos {} -o {} --ref_fasta {} --gc --dinuc_freq {}".format(
            "python",gen_dinucleotide_freqs,basedir + '/test.all.bed',1.0,
            basedir + '/test.final.bed',ref_fasta,basedir + '/test.gc.txt')

    submitter(cmd)
    pos_test = open(basedir + '/test.pos.bed')
    pos_test_len = len(pos_test.readlines())

    os.system('head -n ' + str(pos_test_len) + ' ' + basedir + '/test.final.bed > ' + basedir + '/test.final.pos.bed')
    os.system('tail -n +' + str(pos_test_len + 1) + ' ' + basedir + '/test.final.bed > ' + basedir + '/test.final.neg.bed')

    get_seqs(basedir + '/test.final.pos.bed', basedir + '/test.final.pos.fasta')
    get_seqs(basedir + '/test.final.neg.bed', basedir + '/test.final.neg.fasta')

    pos_test.close()


def main():
    setup_pool()


def setup_pool():
    pool_inputs = []
    clusters=os.listdir(Basedir)
    for cluster in clusters:
        for fold in range(10):
            pool_inputs.append((cluster, fold))
    with ProcessPoolExecutor(max_workers=40) as pool:
        merge=pool.map(train_preprocess, pool_inputs)
    with ProcessPoolExecutor(max_workers=40) as pool:
        merge=pool.map(test_preprocess, pool_inputs)

def run():
    clusters=os.listdir(Basedir)
    for cluster in clusters:
        for fold in range(10):
            stdin=(cluster,fold)
            train_preprocess(stdin)
            test_preprocess(stdin)

if __name__ == "__main__":
    run()
