import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns
import os
import sys
import scipy.stats
import matplotlib as mpl
from itertools import groupby
import math
import scipy.stats
from scipy.stats.mstats import gmean
#from matplotlib_venn import venn2, venn2_circles
from scipy.stats import pearsonr
import random
import itertools


def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()


def NormLog(DF,label,method,O,outDir):
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    def log(L):return np.array([math.log(i+1,2) for i in L])
    if method=='QNorm':
        print('Do QNorm...')
        rank_mean=DF.stack().groupby(DF.rank(method='first').stack().astype(int)).mean()
        DF=DF.rank(method='min').stack().astype(int).map(rank_mean).unstack()
        print('QNorm Done!')
    if method=='DEseq':
        print('Do DEseq...')
        SizeFactor=(DF.T/DF.apply(gmean,axis=1)).dropna(axis=1,how='any').apply(np.median,axis=1)
        DF=DF/SizeFactor
        outSizeFactor=os.path.join(outDir,'DEseq_sizeFactors.txt')
        SizeFactor.to_csv(outSizeFactor,sep='\t')
        print('DEseq Done!')
    DFlog=DF.apply(log)
    if O:
        print('Output file...')
        outfilelog=os.path.join(outDir,'PeakCount.'+label+'_'+method+'_Normalized.log2.txt')
        outfile=os.path.join(outDir,'PeakCount.'+label+'_'+method+'_Normalized.txt')
        DF.to_csv(outfile,sep='\t')
        DFlog.to_csv(outfilelog,sep='\t')
    return DF,DFlog

#Functions
def Mkdir(DirX):
    if not os.path.exists(DirX):
        os.mkdir(DirX)

def Read(File):
    return pd.read_table(File,sep='\t',index_col=0)

def ReadBed(File):
    Bed=pd.read_table(File,sep='\t',index_col=None,header=None)
    Bed.index=Bed[3]
    return Bed

def Save(DF,File):
    DF.to_csv(File,sep='\t')

def SaveTable(DF,FileName):DF.to_csv(FileName,sep='\t')

def SaveBed(Bed,File):
    Bed.to_csv(File,sep='\t',index=False,header=False)
    
def GetState(i):
    if 'HC' in i.upper(): return 'Norm'
    if 'VKH' in i.upper(): return 'VKH'
    if 'GD' in i.upper():return 'GD'
    if 'BD' in i.upper():return 'BD'

def ReadTable(Infile):return pd.read_table(Infile,sep='\t',index_col=0)

def SizeFactor(DF):return (DF.T/DF.apply(gmean,axis=1)).dropna(axis=1,how='any').apply(np.median,axis=1)

def Ttest1(A,B,Index_List):
    return pd.Series(scipy.stats.ttest_ind(A,B,axis=1)[1],index=Index_List)

def Ttest(A,B,Index_List):
    P=[]
    for i in Index_List:
        _,p=scipy.stats.ttest_ind(list(A.loc[i].dropna()),list(B.loc[i].dropna()))
        P.append(p)
    return pd.Series(P,index=Index_List)


def meanCenter(L):
    m=np.mean(L)
    return [i-m for i in L]

def log10(L):return -math.log(L,10)

def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
# input A/B is a DataFrame(log2)
def Diff_FDR(A,B,fdc,fdr):
    Indexs=list(set(list(A[(A==0).apply(sum,axis=1)<len(A.columns)].index)+list(B[(B==0).apply(sum,axis=1)<len(B.columns)].index)))
    A2=A.loc[Indexs]
    B2=B.loc[Indexs]
    MeanA=A2.apply(np.mean,axis=1)
    MeanB=B2.apply(np.mean,axis=1)
    FDC=MeanA-MeanB
    Pval=Ttest(A2,B2,Indexs)
    FDR=pd.Series(p_adjust_bh(Pval.values), index=Indexs)
    return FDC,FDR,list(A2[(FDC>fdc)&(FDR<fdr)].index),list(A2[(FDC<-fdc)&(FDR<fdr)].index)

# input A/B is a DataFrame(log2)
def Diff_FDR(A,B,fdc,fdr):
    Indexs=list(set(list(A[(A==0).apply(sum,axis=1)<len(A.columns)].index)+list(B[(B==0).apply(sum,axis=1)<len(B.columns)].index)))
    A2=A.loc[Indexs]
    B2=B.loc[Indexs]
    MeanA=A2.apply(np.mean,axis=1)
    MeanB=B2.apply(np.mean,axis=1)
    FDC=MeanA-MeanB
    Pval=Ttest(A2,B2,Indexs)
    FDR=pd.Series(p_adjust_bh(Pval.values), index=Indexs)
    return FDC,FDR,list(A2[(FDC>fdc)&(FDR<fdr)].index),list(A2[(FDC<-fdc)&(FDR<fdr)].index)


def Diff_Pval(A,B,fdc,p):
    Indexs=list(set(list(A[(A==0).apply(sum,axis=1)<len(A.columns)].index)+list(B[(B==0).apply(sum,axis=1)<len(B.columns)].index)))
    A2=A.loc[Indexs]
    B2=B.loc[Indexs]
    MeanA=A2.apply(np.mean,axis=1)
    MeanB=B2.apply(np.mean,axis=1)
    FDC=MeanA-MeanB
    Pval=Ttest(A2,B2,Indexs)
    return FDC,Pval,list(A2[(FDC>fdc)&(Pval<p)].index),list(A2[(FDC<-fdc)&(Pval<p)].index)


def GetFDC(CountDF,ct,normal="HC",effect="VKH"):
    NormA=CountDF[[i for i in list(CountDF) if normal in i and ct in i]]
    AffectedA=CountDF[[i for i in list(CountDF) if effect  in i and ct in i]]
    return AffectedA.apply(np.mean,axis=1)-NormA.apply(np.mean,axis=1)

############################
