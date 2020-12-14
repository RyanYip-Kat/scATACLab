import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns
import os
import scipy.stats
from itertools import groupby
import math
import scipy.stats
#from matplotlib_venn import venn2, venn2_circles
import random
import itertools
import sys
import argparse

sys.path.append("/home/ye/Work/BioAligment/SNP/Shi/SSc-ATAC-seq")
from functions import *


############################
print("INFO : get LR pairs")
RLDF1=open('/home/ye/Work/BioAligment/SNP/Shi/SSc-ATAC-seq/SourceData/Receptor-Ligand/DLRP.All_LigandRecepter_Interaction.txt').readlines()
RLDF2=pd.read_table('/home/ye/Work/BioAligment/SNP/Shi/SSc-ATAC-seq/SourceData/Receptor-Ligand/interactions_cellphonedb.csv',sep=',',index_col=0)[['entry_name_a','entry_name_b']].dropna(axis=0,how='any')


A=[i.split('_')[0] for i in RLDF2['entry_name_a']]
B=[i.split('_')[0] for i in RLDF2['entry_name_b']]
RLList=list(set(A+B))
RL=[A[i]+'-'+B[i] for i in range(len(A))]

l=1
for line in RLDF1:
    if line[0]=='#':continue
    if line=='\n':
        l=1
        continue
    if l:
        Gene1=line.split('\t')[1]
        if Gene1 not in RLList:
            RLList.append(Gene1)
        l=0
    else:
        Gene2=line.split('\t')[1]
        if Gene2 not in RLList:
            RLList.append(Gene2)
        if (Gene1+'-'+Gene2 not in RL) and (Gene2+'-'+Gene1 not in RL):
            RL.append(Gene1+'-'+Gene2)

def FoldChange(RL_df,CellType,normal,effect):
        print("INFO : Norm FoldChange")
        #CellType=['CD4','CD8','DC','Fib']
        FDC_df=pd.DataFrame({},index=RL_df.index)
        Header=list(RL_df)
        for ct in CellType:
            Norm=[i for i in Header if (ct in i) and (normal in i)]
            Arm=[i for i in Header if (ct in i) and (effect in i)]
            MNorm=RL_df[Norm].apply(np.mean,axis=1)
            MArm=RL_df[Arm].apply(np.mean,axis=1)
            L=[]
            for i in range(len(MArm)):
                #if max(MArm[i],MNorm[i])<tho:
                #   L.append(0)
                #else:
                L.append(MArm[i]-MNorm[i])
            FDC_df[ct]=L
        return FDC_df

#Find Recepter-Ligand interaction
def FindRLpair(RL,FDC_df_Final):
    RLState=[]
    S=list(FDC_df_Final.index)
    for rl in RL:
        P=rl.split('-')
        if (P[0] in S) and (P[1] in S):
            RLState.append(rl)
    return RLState

def Muxp(a,b):
    if (a<0) and (b<0):
        return min([a,b])
    elif max(a,b)<1:
        return 0
    else:
        return a+b
##Display the Recepter-Ligand（CD4-CD8,CD4-DC,CD4-Fib,CD8-DC,CD8-Fib,DC-Fib）
def GetRLCellTypepair(FDC_df_Final,x,CellType,RLStay):
    RLCellTypes_df=pd.DataFrame({})
    RLCellTypesCount_df=pd.DataFrame({})
    
    #CellType=['CD4','CD8','DC','Fib']
    J=[]
    for i1 in range(len(CellType)-1):
        for i2 in range(i1+1,len(CellType)):
            ct1=CellType[i1]
            ct2=CellType[i2]
            CT1=FDC_df_Final[ct1]
            CT2=FDC_df_Final[ct2]
            S=[]
            C=[]
            RLFinal=[]
            for rl in RLStay:
                try:
                    rl1,rl2=rl.split('-')
                    count=max(Muxp(CT1[rl2],CT2[rl1]),Muxp(CT1[rl1],CT2[rl2]))
                    J.append(count)
                    if count>x:
                        S.append(1)
                        C.append(count)
                    else:
                        S.append(0)
                        C.append(0)
                    RLFinal.append(rl)
                except KeyError:
                    continue
            RLCellTypes_df[ct1+'-'+ct2]=S
            RLCellTypesCount_df[ct1+'-'+ct2]=C
    RLCellTypes_df.index=RLFinal
    RLCellTypesCount_df.index=RLFinal
    RLCellTypes_df=RLCellTypes_df[RLCellTypes_df.apply(sum,axis=1)>1]
    RLCellTypesCount_df=RLCellTypesCount_df.loc[RLCellTypes_df.index]
    return RLCellTypes_df,RLCellTypesCount_df


############################
Rscript="/home/ye/anaconda3/envs/scatac/bin/Rscript"
############################
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Main program')
    parser.add_argument('--count',type=str,default=None,help="scATAC featureCount")
    parser.add_argument('--outdir',type=str,default="./test",help="path to save result")
    parser.add_argument('--annotation',type=str,default=None,help="annotate file from homer")
    parser.add_argument('--normal',type=str,default="HC",help="normal name")
    parser.add_argument('--effect',type=str,default="VKH",help="effect name")
    args = parser.parse_args()

    outdir=args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    name="SSc-scATAC"
    featureCounts=args.count
    outDir=outdir + "/" + name

    normal=args.normal
    effect=args.effect

    print("INFO : Norm Counts")
    CommonDF=Read(featureCounts)
    CommonDF_Norm,CommonDF_Normlog=NormLog(CommonDF,'PBMC','QNorm',True,outDir)  # pass

    CountF=outDir + "/" + "PeakCount.PBMC_QNorm_Normalized.log2.txt"
    CountDF=ReadTable(CountF)
    
    FDC_df=pd.DataFrame({}, index=CountDF.index)
    CellType=list(np.unique([column.split("_")[0] for column in CountDF.columns]))
    for cell in CellType:
       FDC_df[cell]=GetFDC(CountDF,cell,normal,effect)

    OutFDC=outdir + "/" + "FDC_df.txt"
    Save(FDC_df,OutFDC)

    print("INFO : Read homer annotation Peak")
    homer_annotation=args.annotation
    AnoDF=pd.DataFrame(ReadTable(homer_annotation)['Gene Name'])
    
    print("INFO : Get Annotation Counts matrix")
    RLList_st=[]
    Count_df=[]

    for i in RLList:
        if i in list(AnoDF['Gene Name']):
            RLList_st.append(i)
            index=list(AnoDF[AnoDF['Gene Name']==i].index)
            Count_df.append(CountDF.loc[CountDF.index.isin(index)].apply(sum,axis=0))
    Count_df=pd.DataFrame(Count_df,index=RLList_st,columns=list(CountDF))  #  like gene score matrix
    OutCount_df=outdir + "/" + "AnoDF_Count_df.txt"
    Save(Count_df,OutFDC)

    Count_DF=Count_df[Count_df.apply(max,axis=1)>2]
    FDC_df=FoldChange(Count_DF,CellType,normal,effect)
    FDC_df_Final=FDC_df[FDC_df.apply(max,axis=1)>0]
    FDC_df_NormH=(-FDC_df)[(-FDC_df).apply(max,axis=1)>0]

    print("INFO : Find Recepter-Ligand interaction")
    RLStay=FindRLpair(RL,FDC_df_Final)


