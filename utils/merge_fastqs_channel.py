import argparse
import os

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--path",type=str,default="the path of fastqs file")
#parser.add_argument("--r1_txt",type=str,default=None,help="the txt file of R1 fastqs")
#parser.add_argument("--r2_txt",type=str,default=None,help="the txt file of R2 fastqs")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("# The fastqs path is  : {}".format(args.path))
path=args.path
merger="/home/ye/anaconda3/envs/BulkBio/bin/merge_fastq"
fastqs=[file for file in os.listdir(path) if "fastq.gz" in file]

r1=[os.path.join(path,file) for file in fastqs  if "R1" in file]
r2=[os.path.join(path,file) for file in fastqs  if "R2" in file]
assert len(r1)==len(r2)
print("# There are {} pair fastqs".format(len(r1)))

r1_list=" --fastq1 ".join(r1)
r2_list=" --fastq2  ".join(r2)

command="nohup {} --fastq1 {} --fastq2 {} --output-path {} > log 2>&1 & ".format(merger,r1_list,r2_list,args.outdir)
print("# The fastqs merge command is :{}".format(command))
os.system(command)


