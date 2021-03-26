import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd


def get_args():
    parser=argparse.ArgumentParser("Get Embeddings data from SCALE.py Result")
    parser.add_argument("--outdir",type=str,default=None)
    parser.add_argument("--data",type=str,default=None)

    args=parser.parse_args()
    return args

def extractEmbed(adata,outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    assert "latent" in adata.obsm.keys()
    assert "X_umap" in adata.obsm.keys()
    latent = adata.obsm["latent"]
    umap = adata.obsm["X_umap"]
    x_umap = pd.DataFrame(umap,columns=["SCALE#UMAP_1","SCALE#UMAP_2"],index=adata.obs_names.to_list())
    x_latent = pd.DataFrame(latent,columns=["SCALE#latent_"+str(i+1) for i in range(latent.shape[1])],index=adata.obs_names.to_list())
    meta=adata.obs

    umapPath = outdir + "/" + "scale_umap.csv"
    latentPath = outdir + "/" + "scale_latent.csv"
    metaPath = outdir + "/" + "scale_meta.csv"
    print("INFO : UMAP Embed write into  {}\n    Latent Embed write into {}\n   Meta write into {}".format(umapPath,latentPath,metaPath))
    x_latent.to_csv(latentPath,sep=",",index_label="barcode")
    x_umap.to_csv(umapPath,sep=",",index_label="barcode")
    meta.to_csv(metaPath,sep=",",index_label="barcode")


if __name__=="__main__":
    opts=get_args()
    print("********* loading data **********")
    adata=sc.read_h5ad(opts.data)
    outdir=opts.outdir

    print("********* extract embediing *****")
    extractEmbed(adata,outdir)
    
