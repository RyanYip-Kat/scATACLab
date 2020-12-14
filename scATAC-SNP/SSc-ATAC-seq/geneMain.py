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

def Muxp2(a,b):
    if a*b<0:
        return 0
    elif max(abs(a),abs(b))<1:
        return 0
    else:
        return a+b

def GetRLCellTypepair2(FDC_df_Final,x,CellType,AllRL):
    RLCellTypes_df=pd.DataFrame({})
    RLCellTypesCount_df=pd.DataFrame({})
    
    for i1 in range(len(CellType)-1):
        for i2 in range(i1+1,len(CellType)):
            ct1=CellType[i1]
            ct2=CellType[i2]
            CT1=FDC_df_Final[ct1]
            CT2=FDC_df_Final[ct2]
            S=[]
            C=[]
            RLFinal=[]
            for rl in AllRL:
                try:
                    rl1,rl2=rl.split('-')
                    A=Muxp2(CT1[rl2],CT2[rl1])
                    B=Muxp2(CT1[rl1],CT2[rl2])
                    if abs(A)>abs(B):
                        count=A
                    elif abs(B)>abs(A):
                        count=B
                    else:
                        count=0
                    if count>x:
                        S.append(1)
                        C.append(count)
                    elif count<-x:
                        S.append(-1)
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
    RLCellTypes_df=RLCellTypes_df[abs(RLCellTypes_df).apply(sum,axis=1)>1]
    RLCellTypesCount_df=RLCellTypesCount_df.loc[RLCellTypes_df.index]
    return RLCellTypes_df,RLCellTypesCount_df    

def FindRecepterLigand(FDC_df_Final,AllRL,Celltype,outDir,The=7,Thes=1):
    Dir=os.path.join(outDir,'Circus_plot')
    Mkdir(Dir)
    AllIndex=AllRL+[i.split('-')[1]+'-'+i.split('-')[0] for i in AllRL]
    AllCelltype=[Celltype]+[i for i in list(FDC_df_Final) if i!=Celltype]
    NetWorkDF=pd.DataFrame({},index=AllIndex,columns=AllCelltype)
    A=[i.split('-')[0] for i in AllIndex]
    B=[i.split('-')[1] for i in AllIndex]
    NetWorkDF[Celltype]=list(FDC_df_Final.loc[A][Celltype])
    for c in AllCelltype[1:]:
        NetWorkDF[c]=list(FDC_df_Final.loc[B][c])
    def FindFromTo(rl):
        r,l=rl.split('-')
        L=NetWorkDF.loc[rl]
        a=L[0]
        if abs(a)<Thes:
            return [0,0,0,0]
        elif a<0:
            O=[]
            for i in L[1:]:
                if i>-Thes:
                    O.append(0)
                elif a+i<-The:
                    O.append('-'+l)
                else:
                    O.append(0)
            if O==[0,0,0]:
                return [0]+O
            else:
                return ['-'+r]+O
        elif a>0:
            O=[]
            for i in L[1:]:
                if i<Thes:
                    O.append(0)
                elif a+i>The:
                    O.append('+'+l)
                else:
                    O.append(0)
            if O==[0,0,0]:
                return [0]+O
            else:
                return ['+'+r]+O
    for rl in AllIndex:
        L=FindFromTo(rl)
        if L!=[0,0,0,0]:
            NetWorkDF.loc[rl]=L
        else:
            NetWorkDF=NetWorkDF.drop([rl],axis=0)
    SaveTable(NetWorkDF,os.path.join(Dir,Celltype+'_Network.txt'))
    return NetWorkDF,os.path.join(Dir,Celltype+'_Network.txt')


def FindFromTo(FDC_df_Final,AllRL,Celltype,outDir,x=0.2,tho=0.5):
    Index=AllRL+[i.split('-')[1]+'-'+i.split('-')[0] for i in AllRL]
    DF=pd.DataFrame({},index=Index)
    From=[i.split('-')[0] for i in Index]
    To=[i.split('-')[1] for i in Index]
    A=Celltype
    Bs=[i for i in list(FDC_df_Final) if i!=A]
    DF[A]=list(FDC_df_Final.loc[From][A])
    for ct in Bs:
        DF[ct]=list(FDC_df_Final.loc[To][ct])

    for rl in list(DF.index):
        f=rl.split('-')[0]
        t=rl.split('-')[1]
        D=list(DF.loc[rl])
        Judge=[Muxp(D[0],D[1]),Muxp(D[0],D[2]),Muxp(D[0],D[3])]
        if (max(Judge)<x) or (max(D)<tho):
            Dn=[0,0,0,0]
        else:
            if D[0]>0.95:
                Dn=['+'+f]
            else:
                Dn=[f]
            for d in D[1:]:
                if (Muxp(D[0],d)>x):
                    if d>0.95:
                        Dn.append('+'+t)
                    else:
                        Dn.append(t)
                else:
                    Dn.append(0)
        DF.loc[rl]=Dn
    #去除都是0的列
    def judgeZero(L):return [i==0 for i in L]
    #DFjudge=DF.apply(judgeZero,axis=1).apply(sum,axis=1)<4
    DFjudge=DF.apply(judgeZero,axis=1).apply(sum)<4
    Dir=os.path.join(outDir,'Circus_plot_'+str(x))
    Mkdir(Dir)
    Save(DF[DFjudge],os.path.join(Dir,'NetWork_'+Celltype+'.txt')  )
    return DF[DFjudge],os.path.join(Dir,'NetWork_'+Celltype+'.txt')


############################
Rscript="/home/ye/anaconda3/envs/scatac/bin/Rscript"
CircusPlot="/home/ye/Work/BioAligment/SNP/Shi/SSc-ATAC-seq/CircusPlot_ForOneCelltype.R"
############################
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Main program')
    parser.add_argument('--score',type=str,default=None,help="scATAC gene score  featureCount")
    parser.add_argument('--outdir',type=str,default="./test",help="path to save result")
    parser.add_argument('--normal',type=str,default="HC",help="normal name")
    parser.add_argument('--effect',type=str,default="VKH",help="effect name")
    args = parser.parse_args()

    outdir=args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    featureCounts=args.score
    normal=args.normal
    effect=args.effect
    thresh=7

    print("INFO : Read Counts")
    Count_df=Read(featureCounts)
    CellType=list(np.unique([column.split("_")[0] for column in Count_df.columns]))
    
    print("INFO : Get FDC dataframe")
    Count_DF=Count_df[Count_df.apply(max,axis=1)>10]
    FDC_df=FoldChange(Count_DF,CellType,normal,effect)
    FDC_df_Final=FDC_df[FDC_df.apply(max,axis=1)>0]
    FDC_df_NormH=(-FDC_df)[(-FDC_df).apply(max,axis=1)>0]

    print("INFO : Find Recepter-Ligand interaction")
    RLStay=FindRLpair(RL,FDC_df_Final)
    myCmap=sns.light_palette("#CA4641",as_cmap=True)

    print("INFO : Get Recepter-Ligand Celltype Counts matrix")
    RLCellTypes_df,RLCellTypesCount_df=GetRLCellTypepair(FDC_df_Final,thresh,CellType,RLStay)     #more than 7
    Norm_RLCellTypes_df,Norm_RLCellTypesCount_df=GetRLCellTypepair(FDC_df_NormH,thresh,CellType,RLStay)

    print("INFO : Get All Recepter-Ligand")
    AllRL=list(set(list(RLCellTypesCount_df.index)+list(Norm_RLCellTypesCount_df.index)))
    AllRLCellTypes_df,AllRLCellTypesCount_df = GetRLCellTypepair2(FDC_df_Final,thresh,CellType,AllRL)

    SSCup_RowSum=AllRLCellTypes_df[AllRLCellTypes_df>0].fillna(0).apply(sum,axis=1)
    Normup_RowSum=AllRLCellTypes_df[AllRLCellTypes_df<0].fillna(0).apply(sum,axis=1)

    #AllRLCellTypesCount_df=AllRLCellTypesCount_df[['CD8-Fib','CD8-DC','CD4-CD8','CD4-Fib','DC-Fib','CD4-DC']]
    RLOrder=list(AllRLCellTypes_df.apply(sum,axis=1).sort_values().index)
    AllRLCellTypesCount_df=AllRLCellTypesCount_df.loc[RLOrder]

    #SSCup_ColSum=AllRLCellTypes_df[AllRLCellTypes_df>0].fillna(0).apply(sum,axis=0)[['CD8-Fib','CD8-DC','CD4-CD8','CD4-Fib','DC-Fib','CD4-DC']]
    #Normup_ColSum=abs(AllRLCellTypes_df[AllRLCellTypes_df<0].fillna(0).apply(sum,axis=0)[['CD8-Fib','CD8-DC','CD4-CD8','CD4-Fib','DC-Fib','CD4-DC']])
    SSCup_ColSum=AllRLCellTypes_df[AllRLCellTypes_df>0].fillna(0).apply(sum,axis=0)
    Normup_ColSum=abs(AllRLCellTypes_df[AllRLCellTypes_df<0].fillna(0).apply(sum,axis=0))

    print("INFO : Save All Recepter")
    fig1=sns.clustermap(AllRLCellTypesCount_df,figsize=(4,17),cmap='RdBu_r',
            row_cluster=False,col_cluster=False,vmin=-10,vmax=10,yticklabels=AllRLCellTypesCount_df.index,linewidths=0.1,linecolor='w')
    plt.setp(fig1.ax_heatmap.get_yticklabels(), rotation=0,fontsize=10)
    plt.setp(fig1.ax_heatmap.get_xticklabels(), rotation=90,fontsize=20)
    fig1.savefig(os.path.join(outdir,'AllUp_RLCellTypes_df.Cluster.Count.pdf'))
    Save(AllRLCellTypes_df,os.path.join(outdir,'AllUp_RLCellTypes_df.Cluster.Count.txt'))

    SSc_RowSumOrder=SSCup_RowSum.loc[RLOrder]
    Norm_RowSumOrder=Normup_RowSum.loc[RLOrder]
    
    fig2=sns.clustermap(SSc_RowSumOrder,figsize=(0.8,17),cmap='RdBu_r',vmin=-6,vmax=6,
            yticklabels=SSc_RowSumOrder.index,row_cluster=False,col_cluster=False)
    plt.setp(fig2.ax_heatmap.get_yticklabels(), rotation=0,fontsize=10)
    plt.setp(fig2.ax_heatmap.get_xticklabels(), rotation=90,fontsize=20)
    fig2.savefig(os.path.join(outdir,'AllUp.SSc_RLCellTypes_df.sorted_byCount.sumRow.pdf'))
    
    fig3=sns.clustermap(Norm_RowSumOrder,figsize=(0.8,17),cmap='RdBu_r',vmin=-6,vmax=6,
            yticklabels=Norm_RowSumOrder.index,row_cluster=False,col_cluster=False)
    plt.setp(fig3.ax_heatmap.get_yticklabels(), rotation=0,fontsize=10)
    plt.setp(fig3.ax_heatmap.get_xticklabels(), rotation=90,fontsize=20)
    fig3.savefig(os.path.join(outdir,'AllUp.Norm_RLCellTypes_df.sorted_byCount.sumRow.pdf'))
    
    print("INFO : Interaction Network")
    AllRL=list(AllRLCellTypesCount_df.index)
    NetWorkFiles=[]
    for cell in CellType:
        NetWorkDF,NetWorkFile=FindRecepterLigand(FDC_df_Final,AllRL,cell,outdir,The=thresh,Thes=1)
        NetWorkFiles.append(NetWorkFile)
    
    for File in NetWorkFiles:
        print("INFO : Circus Plot for {}".format(File))
        cmd="{} {} --network  {} --palettes {}".format(Rscript,CircusPlot,File,"paired")
        submitter(cmd)

    AllRL=list(RLCellTypes_df.index)
    Dfiles=[]
    for cell in CellType:
        Ddf,Dfile=FindFromTo(FDC_df_Final,AllRL,cell,outdir,10,5)
        Dfiles.append(Dfile)

    AllRL=list(Norm_RLCellTypes_df.index)
    DNormfile=[]
    for cell in CellType:
        Ddf,Dfile=FindFromTo(FDC_df_Final,AllRL,cell,outdir,10,5)
        DNormfile.append(Dfile)



