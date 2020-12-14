import numpy as np
import pandas as pd
import argparse
import os
import subprocess

from sklearn.metrics import roc_auc_score, average_precision_score
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from sklearn.model_selection import StratifiedKFold

####################
rscript="/home/ye/anaconda3/envs/scatac/bin/Rscript"
rcode="/home/ye/Work/BioAligment/SNP/Shi/svm_train/genNullSeqs-LSGKM-One.R"  # --bed --outdir

####################
def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

def makedir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def generate_pos_neg(bed,out):
    cmd="{} {} --bed {} --outdir {}".format(rscript,rcode,bed,out)
    submitter(cmd)

def split_fold(bed_file,nfold=5,output="split"):
    data=pd.read_csv(bed_file,sep="\t",header=None)
    name=os.path.basename(bed_file)

    y=np.random.choice(list(range(0,nfold)),size=data.shape[0],p=[1/nfold]*nfold,replace=True)
    #y=np.random.choice(list(range(0,2)),size=data.shape[0],replace=True)
    X=np.array(data)

    skf = StratifiedKFold(n_splits=nfold)
    i=1

    pos_neg_inputs={}
    for tr_index, ts_index in skf.split(X, y):
        X_train, X_test = X[tr_index], X[ts_index]
        y_train, y_test = y[tr_index], y[ts_index]
        print("INFO ---  Fold : {} --- X_train size : {},X_test size {} ".format(i,X_train.shape[0],X_test.shape[0]))
        X_train=pd.DataFrame(X_train)
        X_test=pd.DataFrame(X_test)

        fold_train_dir=output+"/"+"fold_"+str(i)+"/"+"train"
        fold_test_dir=output+"/"+"fold_"+str(i)+"/"+"test"

        makedir(fold_train_dir)
        makedir(fold_test_dir)
        X_train.to_csv(fold_train_dir+"/"+name,sep='\t', index=False,columns=None,header=False)
        X_test.to_csv(fold_test_dir+"/"+name,sep='\t', index=False,columns=None,header=False)
        i+=1

        inputs=[fold_train_dir+"/"+name,fold_test_dir+"/"+name]
        pos_neg_inputs["fold_"+str(i)]=inputs
        #print("INFO ---  Fold : {} --- generate pos and neg file for train".format(i))
        #cmd="{} {} --bed {} --outdir {}".format(rscript,rcode,fold_train_dir+"/"+name,fold_train_dir)
        #submitter(cmd)

        #print("INFO ---  Fold : {} --- generate pos and neg file for test".format(i))
        #cmd="{} {} --bed {} --outdir {}".format(rscript,rcode,fold_test_dir+"/"+name,fold_test_dir)
        #submitter(cmd)
    return pos_neg_inputs




def get_acc(pos_predfile,neg_predfile,output):
    pos_preds = [float(x.rstrip().split("\t")[1]) for x in open(pos_predfile)]
    neg_preds = [float(x.rstrip().split("\t")[1]) for x in open(neg_predfile)]
    with open(output+'/accuracy.txt', 'w') as acc_file:
        acc_file.write("AUROC: " + str(roc_auc_score(y_true=[1 for x in pos_preds]+[0 for x in neg_preds],
                        y_score = pos_preds+neg_preds)) + '\n')
        acc_file.write("AUPRC: " + str(average_precision_score(y_true=[1 for x in pos_preds]+[0 for x in neg_preds],
                        y_score = pos_preds+neg_preds)) + '\n')
    acc_file.close()


def test_svm(inputs):
    """
    0 : test_seqfile
    1 : model_file
    2 : output_file
    """
    cmd='/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmpredict -T 16 ' + inputs[0] + ' ' + inputs[1] + ' ' + inputs[2]
    submitter(cmd)

def train_svm(inputs):
    """
    0 : posfile
    1 : negfile
    2 : outprefix
    """
    cmd='/home/ye/Work/BioAligment/SNP/lsgkm/bin/gkmtrain -T 16 ' + inputs[0] + ' ' + inputs[1] + ' ' + inputs[2]
    submitter(cmd)


def main(args):
    nfold=args.nfold
    workers=args.workers
    bed_file=args.bed
    output=args.outdir
    split=args.split

    bed_file_dict=split_fold(bed_file,nfold,output)
    train_inputs=[]
    test_pos_inputs=[]
    test_neg_inputs=[]
    for k,v in bed_file_dict.items():
        if not os.path.isfile(os.path.dirname(v[0])+"/"+"ctcfpos.fa"):
            print("INFO --- generate pos and neg file for {}".format(k))
            print("--- generate for train into {}".format(os.path.dirname(v[0])))
            generate_pos_neg(v[0],os.path.dirname(v[0]))

        train_inputs.append((os.path.dirname(v[0])+"/"+"ctcfpos.fa",
            os.path.dirname(v[0])+"/"+"ctcfneg.fa",
            os.path.dirname(v[0])+"/"+"train"))

        if not os.path.isfile(os.path.dirname(v[1])+"/"+"ctcfpos.fa"):
            print("--- generate for test into {}".format(os.path.dirname(v[1])))
            generate_pos_neg(v[1],os.path.dirname(v[1]))

        test_pos_inputs.append((os.path.dirname(v[1])+"/"+"ctcfpos.fa",
            os.path.dirname(v[0])+"/"+"train.model.txt",
            os.path.dirname(v[1])+"/"+"test.pos.pred"))

        test_neg_inputs.append((os.path.dirname(v[1])+"/"+"ctcfneg.fa",
            os.path.dirname(v[0])+"/"+"train.model.txt",
            os.path.dirname(v[1])+"/"+"test.neg.pred"))
    
    print("INFO : TRAIN")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        merge=pool.map(train_svm,train_inputs)
    
    print("INFO : PREDICT POS")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        merge=pool.map(test_svm,test_pos_inputs)
    
    print("INFO : PREDICT NEG")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        merge=pool.map(test_svm,test_neg_inputs)

    print("INFO : get accuracy")
    for k,v in bed_file_dict.items():
        try:
            pos=os.path.dirname(v[1])+"/"+"test.pos.pred"
            neg=os.path.dirname(v[1])+"/"+"test.neg.pred"
            out=os.path.dirname(v[1])
            print("INFO --- Pos : {},Neg {},Out : {}".format(pos,neg,out))
            get_acc(pos,neg,out)
        except:
            print("Please check your pos and neg file!!!")


    #return train_inputs,test_pos_inputs,test_neg_inputs

    
def get_args():
    parser = argparse.ArgumentParser("CV train and test")
    parser.add_argument('--nfold',type=int,default=5, help='the number fold')
    parser.add_argument('--workers',type=int,default=12, help='the number workers')
    parser.add_argument('--bed',type=str,default=None, help='the bed file use for training and testing model')
    parser.add_argument('--outdir',type=str,default="./split", help='the path to save result')
    parser.add_argument('--split',action="store_true", help='wether to split dataset')
    args = parser.parse_args()
    return args

if __name__=="__main__":
    args=get_args()
    main(args)

