library(argparse)
library(stringr)
library(data.table)

#############################
parser <- ArgumentParser(description='Process some tasks')
parser$add_argument("--fragment",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="fragment files")


args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

bgzip="/home/ye/anaconda3/envs/BulkBio/bin/bgzip"
tabix="/home/ye/anaconda3/envs/BulkBio/bin/tabix"

#############################
ZipAndBuildTbi=function(gzipFile){
        fileName=file.path(dirname(gzipFile),paste0("sorted-",basename(gzipFile)))
        cmd=paste0("sort -k1,1d -k2,2n -k3,3n"," ",gzipFile," > ",fileName)
        system(cmd)
        cmd=paste0(bgzip," ",fileName)
        message("INFO : Zip File ...")
        system(cmd)

        cmd=paste0(tabix," ","--preset=bed"," ",fileName,".gz")
        message("INFO : Build tabix Index  ...")
        system(cmd)
        system(paste0("rm ",gzipFile))
        return(paste0(fileName,".gz"))
}


fragment_list=args$fragment
for(i in seq_along(fragment_list)){
	fragFile=fragment_list[i]
	cat(sprintf("INFO : Loading No. [ %d of %d ] fragmnet file from --- [ %s ]\n",i,length(fragment_list),fragFile))
	DF=fread(fragFile)
	message("INFO : Make new barcode...")
	DF$V4=paste(DF$V4,1:nrow(DF),sep="-")
	fragOut=file.path(dirname(fragFile),paste0("reformat-",basename(fragFile)))
	message("INFO : Write ...")
	write.table(DF,fragOut,sep="\t",quote=F,row.names=F,col.names=F)
	ZipAndBuildTbi(fragOut)
}

message("INFO : Done!")

