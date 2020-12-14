import argparse
import subprocess
import numpy as np
import os
import vizsequence
import modisco
import pickle
import h5py

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
    parser.add_argument('--fasta',type=str, help='the squence fasta file (which can be genrated with bedtools getfasta)',required=True)
    parser.add_argument('--model',type=str, help='the model file from gkmtrain',required=True)
    parser.add_argument('--outdir',type=str,help='the output path',default="./results")
    parser.add_argument('--gkmexplain',type=str,help='the path of gkmexplain',default='/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmexplain')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    print("INFO : Run gkmexplain")
    explain_result=os.path.join(args.outdir,"explanation.txt")
    explain_cmd="{} {} {} {} {}".format(args.gkmexplain,"-m 0",args.fasta,args.model,explain_result)
    submitter(explain_cmd)

    hyp_explain_result=os.path.join(args.outdir,"hyp_explanation.txt")
    hyp_explain_cmd="{} {} {} {} {}".format(args.gkmexplain,"-m 1",args.fasta,args.model,hyp_explain_result)
    submitter(hyp_explain_cmd)


    print("INFO : Get score")
    impscores = [
            np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
            for x in open(explain_result)]
    
    hyp_impscores = [
            np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
            for x in open(hyp_explain_result)]

    fasta_seqs = [x.rstrip() for (i,x) in enumerate(open(args.fasta))
              if i%2==1]
    
    onehot_data = np.array([one_hot_encode_along_channel_axis(x)
                         for x in fasta_seqs])

    for i in range(len(impscores)):
        #viz_sequence.plot_weights(impscores[i], subticks_frequency=10)
        plot_weights(impscores[i], subticks_frequency=10,ylabel="sequence_"+str(i))
        filename=os.path.join(args.outdir,"sequence_"+str(i)+"_impscores.pdf")
        plt.savefig(filename)

        imp_score_each_pos = np.sum(impscores[i],axis=-1)
        imp_score_sign_each_pos = np.sign(imp_score_each_pos)
        hyp_scores_same_sign_mask = (np.sign(hyp_impscores[i])*imp_score_sign_each_pos[:,None] > 0)

        hyp_scores_same_sign_imp_scores_sum = np.sum(hyp_impscores[i]*hyp_scores_same_sign_mask,axis=-1)
        norm_ratio = imp_score_each_pos/hyp_scores_same_sign_imp_scores_sum
        norm_hyp = hyp_impscores[i]*norm_ratio[:,None]
        
        plot_weights(norm_hyp, subticks_frequency=10)
        filename=os.path.join(args.outdir,"sequence_"+str(i)+"_norm_hyp.pdf")
        plt.savefig(filename)

        plot_weights(norm_hyp*onehot_data[i], subticks_frequency=10)
        filename=os.path.join(args.outdir,"sequence_"+str(i)+"_norm_hyp_onehot.pdf")
        plt.savefig(filename)

    normed_hyp_scores = []
    normed_impscores = []
    for i in range(len(impscores)):
        imp_score_each_pos = np.sum(impscores[i],axis=-1)
        imp_score_sign_each_pos = np.sign(imp_score_each_pos)
        hyp_scores_same_sign_mask = (np.sign(hyp_impscores[i])*imp_score_sign_each_pos[:,None] > 0)
        hyp_scores_same_sign_imp_scores_sum = np.sum(hyp_impscores[i]*hyp_scores_same_sign_mask,axis=-1)
        norm_ratio = imp_score_each_pos/hyp_scores_same_sign_imp_scores_sum
        norm_hyp = hyp_impscores[i]*norm_ratio[:,None]
        normed_hyp_scores.append(norm_hyp)
        normed_impscores.append(norm_hyp*onehot_data[i])


    tfmodisco_results = modisco.tfmodisco_workflow.workflow.TfModiscoWorkflow(
            sliding_window_size=6,
            flank_size=2,
            seqlets_to_patterns_factory=modisco.tfmodisco_workflow.seqlets_to_patterns.TfModiscoSeqletsToPatternsFactory(
                initclusterer_factory=modisco.clusterinit.memeinit.MemeInitClustererFactory(
                        meme_command="meme", base_outdir="meme_out",
                        #max_num_seqlets_to_use specifies the maximum number of seqlets to use
                        # with MEME (this is to speed up MEME in the cases where the number of seqlets is
                        #  very large)
                        max_num_seqlets_to_use=10000,
                        nmotifs=10,
                        n_jobs=4),
                trim_to_window_size=6,
                initial_flank_to_add=2,
                final_flank_to_add=5,
                kmer_len=6, num_gaps=1,
                num_mismatches=0))(
                task_names=["task0"],
                contrib_scores={'task0': normed_impscores},                
                hypothetical_contribs={'task0': normed_hyp_scores},
                one_hot=onehot_data)

    #filename=open(os.path.join(args.outdir,"tfmodisco_results.pkl"),"w")
    #pickle.dump(tfmodisco_results,filename)

    grp = h5py.File(os.path.join(args.outdir,"tfmodisco_results.hdf5"),"w")
    tfmodisco_results.save_hdf5(grp)
