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


parser$add_argument("--batch_correct",
		    action="store_true",
		    default=FALSE,
		    help="whether run batch correct")

parser$add_argument("--column",
                    type="character",
                    default=NULL,
		    help="column in snap metaData")

parser$add_argument("--subset",
		    nargs="+",
		    type="character",
                    default=NULL)
args <- parser$parse_args()

#############################
outDir = args$outdir
Temp=file.path(outDir,"tmp")
mkdir(outDir)
mkdir(Temp)


#############################
message("INFO : loading dataset ...")
x.sp=readRDS(args$snap)
metadata=x.sp@metaData
if(!is.null(args$column) & !is.null(args$subset)){
	assertthat::assert_that(args$column%in%colnames(metadata))
	metadata[[args$column]]=as.character(metadata[[args$column]])
	column=metadata[[args$column]]
	idx=args$subset%in%unique(column)
	assertthat::assert_that(sum(idx)==length(args$subset))
	idx=which(column%in%args$subset)
	x.sp=x.sp[idx,]
}

############################# dataset size < 20k
message("INFO : run DiffusionMap ...")
if(length(x.sp@barcode)<20000){
       x.sp = runDiffusionMaps(
			       obj=x.sp,
			       input.mat="bmat",
			       num.eigs=50);
}else{
	############################### datatset > 20k,use this projection method
        row.covs = log10(Matrix::rowSums(x.sp@bmat)+1);
        row.covs.dens = density(
				x = row.covs,
				bw = 'nrd', adjust = 1
				);

        sampling_prob = 1 / (approx(x = row.covs.dens$x, y = row.covs.dens$y, xout = row.covs)$y + .Machine$double.eps);
        set.seed(1);
        idx.landmark.ds = sort(sample(x = seq(nrow(x.sp)), size = 15000, prob = sampling_prob));
        x.landmark.sp = x.sp[idx.landmark.ds,];
        x.query.sp = x.sp[-idx.landmark.ds,];
	x.landmark.sp = runDiffusionMaps(
					 obj= x.landmark.sp,
					 input.mat="bmat",
					 num.eigs=50);
        x.query.sp = runDiffusionMapsExtension(
					       obj1=x.landmark.sp,
					       obj2=x.query.sp,
					       input.mat="bmat");
        x.landmark.sp@metaData$landmark = 1;
        x.query.sp@metaData$landmark = 0;
        x.sp = snapRbind(x.landmark.sp, x.query.sp);
        ## combine landmarks and query cells;
        x.sp = x.sp[order(x.sp@sample),]; # IMPORTANT
        rm(x.landmark.sp, x.query.sp); # free memory
}

plotDimReductPW(obj=x.sp, 
    eigs.dims=1:50,
    point.size=0.3,
    point.color="grey",
    point.shape=19,
    point.alpha=0.6,
    down.sample=5000,
    pdf.file.name=file.path(outDir,"plotDimReductPW.pdf"), 
    pdf.height=10, 
    pdf.width=12
  );

############################# remove batch correct
if(args$batch_correct){
	message("INFO : remove batch correct ...")
	x.sp = runHarmony(
				obj=x.sp,
				eigs.dim=1:20,
				meta_data=x.sp@sample # sample index
				);
}

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
     resolution=1.2,
     seed.use=10
  );
x.sp@metaData$cluster = x.sp@cluster;

message("INFO : run reduction ...")
message("INFO : run tsne ...")
x.sp = runViz(
    obj=x.sp,
    tmp.folder=Temp,
    dims=2,
    eigs.dims=1:20,
    method="Rtsne",
    seed.use=10,
    num.cores=16
  );

message("INFO : run umap ...")
x.sp = runViz(
    obj=x.sp,
    tmp.folder=Temp,
    dims=2,
    eigs.dims=1:20,
    method="umap",
    seed.use=10,
    num.cores=16,
  );


message("INFO : Save ...")
saveRDS(x.sp,file.path(outDir,"snapCluster.rds"))
message("INFO : Done!")
