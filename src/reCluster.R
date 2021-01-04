library(SnapATAC)
library(harmony)
library(GenomicRanges)
library(argparse)
library(stringr)
############################
mkdir=function(path){
	if(!dir.exists(path)){
		dir.create(path,recursive=TRUE)
	}
}

findMacs2=function(){
        macs2_path=system("which macs2",intern=TRUE)
        return(macs2_path)
}

findSnaptools=function(){
        path=system("which snaptools",intern=TRUE)
        return(path)
}

############################
parser <- ArgumentParser(description='Create snap Object from mutiple sample ...')
parser$add_argument("--snap",
                    type="character",
                    default=NULL,
                    help="snap object file from create Object")



parser$add_argument("--outdir",
                    type="character",
                    default="Snap_result")


parser$add_argument("--resolution",
		    type="double",
                    default=1.2,
		    help="resolution value for cluster")
args <- parser$parse_args()

#############################
outDir = args$outdir
Temp=file.path(outDir,"tmp")
resolution=args$resolution
mkdir(outDir)
mkdir(Temp)


#############################
message("INFO : loading dataset ...")
x.sp=readRDS(args$snap)

############################## Graph-based cluster
message("INFO : Graph-based cluster ...")
x.sp = runKNN(
    obj= x.sp,
    eigs.dim=1:20,
    k=15
  );
x.sp = runCluster(
     obj=x.sp,
     tmp.folder=Temp,
     louvain.lib="R-igraph",
     path.to.snaptools=findSnaptools(),
     resolution=resolution,
     seed.use=10
  );
newColumn=paste0("cluster_",resolution)
x.sp@metaData[[newColumn]] = x.sp@cluster;
message("INFO : Save ...")
saveRDS(x.sp,file.path(outDir,"snapCluster.rds"))
message("INFO : Done!")
