import argparse
import subprocess
import numpy as np
import os
import vizsequence

from vizsequence import viz_sequence
from matplotlib import pyplot as plt

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

def one_hot_encode_along_channel_axis(sequence):
    to_return = np.zeros((len(sequence),4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return,
                                 sequence=sequence, one_hot_axis=1)
    return to_return

def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    assert one_hot_axis==0 or one_hot_axis==1
    if (one_hot_axis==0):
        assert zeros_array.shape[1] == len(sequence)
    elif (one_hot_axis==1): 
        assert zeros_array.shape[0] == len(sequence)
    #will mutate zeros_array
    for (i,char) in enumerate(sequence):
        if (char=="A" or char=="a"):
            char_idx = 0
        elif (char=="C" or char=="c"):
            char_idx = 1
        elif (char=="G" or char=="g"):
            char_idx = 2
        elif (char=="T" or char=="t"):
            char_idx = 3
        elif (char=="N" or char=="n"):
            continue #leave that pos as all 0's
        else:
            raise RuntimeError("Unsupported character: "+str(char))
        if (one_hot_axis==0):
            zeros_array[char_idx,i] = 1
        elif (one_hot_axis==1):
            zeros_array[i,char_idx] = 1

def plot_weights(array,
                 figsize=(20,2),
                 ax_transform_func=lambda x: x,
                 **kwargs):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111) 
    ax_transform_func(ax)
    viz_sequence.plot_weights_given_ax(ax=ax,
        array=array,
        **kwargs)
    
def deltasvm_scores(sequence, deltasvm_kmer_to_pred, lmersize):
  
  onehot_sequence = np.zeros((len(sequence),4))
  scores = np.zeros((len(sequence),4))
  for (i,base) in enumerate(sequence):
    if (i <= len(sequence)-lmersize):
      lmer = sequence[i:i+lmersize]
      lmer_pred = deltasvm_kmer_to_pred[lmer]
      for (j,lmer_base) in enumerate(lmer):
        for base_idx,base_sub in enumerate(['A','C','G','T']):
          if base_sub!=lmer_base:
            new_lmer = lmer[:j]+base_sub+lmer[(j+1):]
            new_pred = deltasvm_kmer_to_pred[new_lmer]
            delta = new_pred - lmer_pred
            scores[i+j, base_idx] = delta
          else:
            onehot_sequence[i+j, base_idx] = 1
  scores = scores - np.mean(scores,axis=1)[:,None]
  return scores*onehot_sequence, scores


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='A program to visualize gkm svm explain result')
    parser.add_argument('--kmer',type=str, help='the kmer squence fasta generated from nrkmers.py')
    parser.add_argument('--model',type=str, help='the model file from gkmtrain',required=True)
    parser.add_argument('--outdir',type=str,help='the output path',default="./results")
    parser.add_argument('--gkmpredict',type=str,help='the path of gkmexplain',default='/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmpredict')
    parser.add_argument('--sequence',type=str,help='the sequence')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    print("INFO : Run gkmpredict")
    explain_result=os.path.join(args.outdir,"preds_kmers.txt")
    explain_cmd="{} {} {} {}".format(args.gkmpredict,args.kmer,args.model,explain_result)
    submitter(explain_cmd)
    
    print("INFO : Run deltasvm_scores")
    deltasvm_kmer_to_pred = {}
    rc_trans = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    for line in open(os.path.join(args.outdir,"preds_kmers.txt")):
        kmer,pred = line.rstrip().split("\t")
        pred = float(pred)
        deltasvm_kmer_to_pred[kmer] = pred
        deltasvm_kmer_to_pred["".join([rc_trans[x] for x in kmer[::-1]])] = pred

    the_seq=args.sequence
    delta_svm_scores, delta_svm_hypscores = deltasvm_scores(the_seq, deltasvm_kmer_to_pred, 6)
    
    plot_weights(delta_svm_scores, subticks_frequency=10)
    filename=os.path.join(args.outdir,"delta_svm_scores.pdf")
    plt.savefig(filename)

    plot_weights(delta_svm_hypscores, subticks_frequency=10)
    filename=os.path.join(args.outdir,"delta_svm_hypscores.pdf")
    plt.savefig(filename)


