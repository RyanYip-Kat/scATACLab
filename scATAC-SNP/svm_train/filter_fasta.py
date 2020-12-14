import os
import argparse
import pickle

chroms=["chr"+str(i) for i in range(1,23)]+["chrX","chrY","chrM","chrx","chry","chrm"]
def skip_fasta(fasta,pattern="N",outdir="./outdir"):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fi=open(fasta)
    lines=fi.readlines()
    fi.close()

    name_idx=[ i for i,line in enumerate(lines) if ">" in line.rstrip()]
    name_idx2=name_idx[1:]+[len(lines)]
    
    assert len(name_idx)==len(name_idx2)
    
    newfasta=outdir+"/"+os.path.basename(fasta)
    fn=open(newfasta,"w")
    seq_dict={}
    for start,end in zip(name_idx,name_idx2):
        the_seqs=[line.rstrip() for line in lines[start+1:end]]
        the_seq="".join(the_seqs)
        name=lines[start]
        chrom=name.split(":")[0].split(">")[1]
        if pattern not in the_seq and chrom in chroms:
            #print("INFO : write name {}".format(name.rstrip()))
            fn.write(name)
            for line in lines[start+1:end]:
                fn.write(line)

            seq_dict[name.rstrip()]=the_seq
    fn.close()

    newDict=outdir+"/"+"sequence_dict.pkl"
    fn=open(newDict,"wb")
    pickle.dump(seq_dict,fn)
    fn.close()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Program to filter fasta file")
    parser.add_argument("fasta")
    parser.add_argument("outdir")
    parser.add_argument("--pattern",type=str,default="N",help="the patther to detect in the sequence")
    args=parser.parse_args()

    skip_fasta(args.fasta,args.pattern,args.outdir)
