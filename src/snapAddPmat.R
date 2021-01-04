library(argparse)
library(SnapATAC)
library(stringr)

mkdir=function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}



##############################
parser <- ArgumentParser(description='Create snap Object from mutiple sample ...')
parser$add_argument("--snap",
                    type="character",
                    default=NULL,
                    help="snap object file from create Object")



args <- parser$parse_args()

#############################


############################
message("INFO : loading dataset ...")
x.sp=readRDS(args$snap)

message("INFO : addPmat ...")
x.sp=addPmatToSnap(x.sp,do.par=TRUE,num.cores=16)

message("INFO : Save ...")
saveRDS(x.sp,args$snap)
message("INFO : Done!")
