import os
import argparse
import numpy as np


def run_ism(model_file_path,sequence,temp_filepath):
  
  letter_to_idx = {'A':0, 'C':1, 'G':2, 'T':3, 'N': -1}
  
  fh = open(temp_filepath,'w')
  fh.write(">orig_seq\n")
  fh.write(sequence+"\n")
  fh.close()

  explanation_file = temp_filepath+".explanation.txt"
  os.system("/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmexplain "+temp_filepath+" "+model_file_path+" "+explanation_file)
  gkmexplain_impscores = [
    np.array( [[float(z) for z in y.split(",")] for y in x.rstrip().split("\t")[2].split(";")])
    for x in open(explanation_file)
  ][0]
  os.remove(temp_filepath)

  seq_len = len(sequence)
  fh = open(temp_filepath,'w')
  fh.write(">orig_seq\n")
  fh.write(sequence+"\n")

  onehot_seq = np.zeros((seq_len,4))  
  for pos in range(seq_len):
    orig_letter_idx = letter_to_idx[sequence[pos]]
    if (orig_letter_idx != -1):
      onehot_seq[pos,orig_letter_idx] = 1
    for letter in ['A','C','G','T']:
      #only need to compute scores for perturbations
      if (letter != sequence[pos]):
        letter_idx = letter_to_idx[letter]
        fh.write(">pos-"+str(pos)+"_base-"+str(letter_idx)+"\n")
        fh.write(sequence[:pos]+letter+sequence[pos+1:]+"\n")
  fh.close()
  
  predictions_file = temp_filepath+".preds.txt"
  os.system("/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmpredict -v 1 "+temp_filepath+" "+model_file_path+" "+predictions_file)

    
  ism_importance_scores = np.zeros((seq_len,4))
  for line in open(predictions_file):
    seq_name,pred = line.rstrip().split("\t")
    pred = float(pred)
    if (seq_name=="orig_seq"):
      orig_pred = pred
    else:
      pos,letter_idx = seq_name.split("_")
      pos = int(pos.split("-")[1])
      letter_idx = int(letter_idx.split("-")[1])
      ism_importance_scores[pos,letter_idx] = (pred - orig_pred)
  
  os.remove(temp_filepath)
  os.remove(predictions_file)
  
  ism_hyp_importance_scores = ism_importance_scores-np.mean(ism_importance_scores,axis=1)[:,None]
  ism_importance_scores = ism_hyp_importance_scores*onehot_seq
  
  return ism_importance_scores, ism_hyp_importance_scores,gkmexplain_impscores


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="the program to run ism score")
    parser.add_argument("model_file_path")
    parser.add_argument("input_file")
    parser.add_argument("outdir")
    args=parser.parse_args()

    model_file_path=args.model_file_path   # "params_t3_l6_k5_d1_g2_c10_w3.model.txt"
    input_file=args.input_file    # "test_positives.fa"
    outdir=args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ism_imp_scores = []
    for i,x in enumerate(open(input_file)):
        if (i%2==1):
            ism_imp_scores.append(
            run_ism(model_file_path=model_file_path,
                    sequence=x.rstrip(),
                    temp_filepath="ismtmpfile.txt")[0])
        if ((i-1)%10==0):
            print("Done",int((i-1)/2))

    np.save(outdir+"/"+"ism_imp_scores", np.array(ism_imp_scores))
    # plot_weights(ism_scores, subticks_frequency=10)





