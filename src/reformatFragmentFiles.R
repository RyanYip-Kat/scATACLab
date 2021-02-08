library(argparse)
library(stringr)
library(Seurat)
library(ArchR)
library(Matrix)
library(ggplot2)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--inputFiles",
                    type="character",
                    default="")


parser$add_argument("--genome",
                    type="character",
                    default="hg38",
		    choices=c("hg38","hg19","mm9","mm10"),
		    help="which genome to be used")

parser$add_argument("--num_threads",
                    type="integer",
                    default=16,
                    help="number of threads for command")

args <- parser$parse_args()


inputFiles<-read.table(args$inputFiles,sep=",",stringsAsFactors=FALSE)
files=as.character(inputFiles$V1)
filenames=as.character(inputFiles$V2)

print(paste0("Setting default genome  with :",args$genome))
addArchRGenome(args$genome)

print(paste0("Setting threads  :",args$num_threads))
addArchRThreads(threads = args$num_threads) 

print("# reformatFragmentFiles")

for(file in files){
	print(paste0("##  Reformat FragmenetsFiles : ",file))
	reformatFragmentFiles(file)
	#format_file=file.path(dirname(file),"fragments-Reformat.tsv.gz")
	#commad=paste0("/home/ye/anaconda3/envs/scATAC/bin/tabix"," ","-p bed -f"," ",format_file)
	#print(commad)
        #system(commad)
}	
