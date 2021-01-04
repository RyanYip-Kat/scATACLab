library(SnapATAC)
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



args <- parser$parse_args()

############################
outDir=args$outdir
mkdir(outDir)

###########################
message("INFO : loading dataset ...")
x.sp=readRDS(args$snap)

########################### plot
plotViz(
    obj=x.sp,
    method="tsne",
    main="Cluster",
    point.color=x.sp@cluster,
    point.size=0.1,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    legend.add=TRUE,
    pdf.file.name=file.path(outDir,"tsne-cluster.pdf"),
    pdf.width=12,
    pdf.height=10
  );

plotViz(
    obj=x.sp,
    method="tsne",
    main="Cluster",
    point.color=x.sp@sample,
    point.size=0.1,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    legend.add=TRUE,
    pdf.file.name=file.path(outDir,"tsne-sample.pdf"),
    pdf.width=12,
    pdf.height=10
  );


plotViz(
    obj=x.sp,
    method="umap",
    main="Cluster",
    point.color=x.sp@cluster,
    point.size=0.1,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    legend.add=TRUE,
    pdf.file.name=file.path(outDir,"umap-cluster.pdf"),
    pdf.width=12,
    pdf.height=10
  );

plotViz(
    obj=x.sp,
    method="umap",
    main="Cluster",
    point.color=x.sp@sample,
    point.size=0.1,
    text.add=TRUE,
    text.size=1,
    text.color="black",
    text.halo.add=TRUE,
    text.halo.color="white",
    text.halo.width=0.2,
    down.sample=10000,
    legend.add=TRUE,
    pdf.file.name=file.path(outDir,"umap-sample.pdf"),
    pdf.width=12,
    pdf.height=10
  );

