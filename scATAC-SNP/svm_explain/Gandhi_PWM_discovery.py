import argparse
import os
import numpy as np

def run_me(model_file_path,nrkmers_fasta,outdir):
    nrkmers_scores=outdir+"/"+"nrkmers_scores.txt"
    cmd="/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmpredict"+" " + "-v 1" + \
             " " + nrkmers_fasta + " " + model_file_path + " " + nrkmers_scores
    print("INFO : Predict kmers fasta")
    os.system(cmd)

    print("INFO : Run gandhietalpwms")
    gandhietalpwms=outdir+"/"+"gandhietalpwms"
    cmd="/home/ye/anaconda3/envs/ldsc/bin/python" + " " + "/home/ye/Work/BioAligment/SNP/gkmsvm-2.0/scripts/svmw_emalign.py" + \
            " " + nrkmers_scores + " " +  str(19) + " " + gandhietalpwms
    os.system(cmd)
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Program to run : PWM discovery method in Gandhi et al. (2014)")
    parser.add_argument("model_file_path")
    parser.add_argument("nrkmers_fasta")
    parser.add_argument("outdir")
    args=parser.parse_args()

    model_file_path=args.model_file_path
    nrkmers_fasta=args.nrkmers_fasta
    outdir=args.outdir

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    run_me(model_file_path,nrkmers_fasta,outdir)



